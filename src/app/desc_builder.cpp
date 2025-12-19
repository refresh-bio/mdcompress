#include "desc_builder.h"
#include <chemfiles.hpp>
#include <queue>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::string;
using std::queue;

//for internal linkage
namespace {
    // Simple hash for (name, natom) key
    struct MolKey {
        std::string name;
        size_t natom;
        bool operator==(const MolKey& other) const {
            return natom == other.natom && name == other.name;
        }
    };

    struct MolKeyHash {
        std::size_t operator()(const MolKey& k) const {
            std::hash<std::string> hs;
            std::hash<size_t>      hi;
            return hs(k.name) ^ (hi(k.natom) + 0x9e3779b97f4a7c15ULL + (hs(k.name) << 6) + (hs(k.name) >> 2));
        }
    };

    struct internal_segment_desc
    {
	    enum class segment_type {mol, wat, other, lip};
        segment_type type;
        uint32_t size;
        int group_id{};
    };

    class Classifier
    {
        std::unordered_map<std::string, int> mol_map;
        std::unordered_map<std::string, int> wat_map;
        std::unordered_map<std::string, int> other_map;
        std::unordered_map<std::string, int> lip_map;

        int get_group_id(std::unordered_map<std::string, int>& map, const std::string& mol_name)
        {
            auto it = map.find(mol_name);
            if (it != map.end())
                return it->second;
            return map.insert(std::make_pair(mol_name, map.size())).first->second;
        }

    public:
        internal_segment_desc classify_group(size_t mol_count,
            size_t natom,
            const std::string& mol_name) {
            internal_segment_desc res;
            res.size = static_cast<uint32_t>(natom);

            if (natom == 3 && mol_count > 100) {
                res.type = internal_segment_desc::segment_type::wat;
                res.group_id = get_group_id(wat_map, mol_name);
            }
            else if (natom == 1) {
                res.type = internal_segment_desc::segment_type::other;
                res.group_id = get_group_id(other_map, mol_name);
            }
            else if (natom > 1000) {
                res.type = internal_segment_desc::segment_type::mol;
                res.group_id = get_group_id(mol_map, mol_name);
            }
            else if (natom > 30 && natom < 120 && mol_count > 50) {
                res.type = internal_segment_desc::segment_type::lip;
                res.group_id = get_group_id(lip_map, mol_name);
            }
            else if (natom > 3 && natom < 300) {
                res.type = internal_segment_desc::segment_type::mol;
                res.group_id = get_group_id(mol_map, mol_name);
            }
            else {
                throw std::runtime_error(
                    "unknown type when classifying atoms, #atoms = " +
                    std::to_string(natom) + ", mol_count = " + std::to_string(mol_count)
                );
            }

            return res;
        }
    };

    // Pick a "molecule name" from residues/atoms
    std::string guess_molecule_name(const chemfiles::Topology& topo,
        const std::vector<size_t>& atoms_in_mol) {
        // 1. Try residue names (e.g. SOL / LIG / PROA)
        std::unordered_map<std::string, size_t> res_name_count;

        for (auto idx : atoms_in_mol) {
            auto opt_res = topo.residue_for_atom(idx);
            if (opt_res) {
                const chemfiles::Residue& res = *opt_res;
                auto name = res.name();
                if (!name.empty()) {
                    res_name_count[name]++;
                }
            }
        }
        if (!res_name_count.empty()) {
            // Choose most frequent residue name in this molecule
            auto best = std::max_element(
                res_name_count.begin(), res_name_count.end(),
                [](auto& a, auto& b) { return a.second < b.second; }
            );
            return best->first;
        }

        // 2. Fallback: first atom type or element symbol
        if (!atoms_in_mol.empty()) {
            const auto& atom = topo[atoms_in_mol.front()];
            if (!atom.name().empty()) {
                return atom.name();
            }
            if (!atom.type().empty()) {
                return atom.type();
            }
        }

        return "UNK";
    }
}

vector<segment_desc_t> build_desc(const string& trajectory_path)
{
    try {
        chemfiles::Trajectory tpr(trajectory_path);
        auto frame = tpr.read();
        const auto& topo = frame.topology();
        const size_t natoms = topo.size();
        if (natoms == 0) {
            throw std::runtime_error("No atoms in trajectory file - " + trajectory_path);
        }
        if (topo.bonds().empty()) {
            throw std::runtime_error("No bonds in topology – " + trajectory_path +
                " is probably not a topology file (TPR/PSF/etc.)");
        }

        // --- Build adjacency list from bonds (to find molecules as connected components)
        std::vector<std::vector<size_t>> adj(natoms);
        for (const auto& bond : topo.bonds()) {
            auto i = bond[0];
            auto j = bond[1];
            if (i >= natoms || j >= natoms)
                throw std::runtime_error("Error while building desc: i >= natoms || j >= natoms, for i = "
                    + std::to_string(i)
                    + ", j = "
                    + std::to_string(j)
                    + ", natoms = "
                    + std::to_string(natoms)
                );
            adj[i].push_back(j);
            adj[j].push_back(i);
        }

        // --- Connected components = molecules
        std::vector<int> comp_id(natoms, -1);
        std::vector<std::vector<size_t>> components;
        int current_comp = 0;

        for (size_t start = 0; start < natoms; ++start) {
            if (comp_id[start] != -1) continue;

            // BFS
            std::queue<size_t> q;
            q.push(start);
            comp_id[start] = current_comp;

            std::vector<size_t> atoms_in_comp;
            atoms_in_comp.push_back(start);

            while (!q.empty()) {
                size_t u = q.front();
                q.pop();
                for (auto v : adj[u]) {
                    if (comp_id[v] == -1) {
                        comp_id[v] = current_comp;
                        q.push(v);
                        atoms_in_comp.push_back(v);
                    }
                }
            }

            components.push_back(std::move(atoms_in_comp));
            current_comp++;
        }

        // --- Group molecules into molblocks by (guessed_name, natom)
        std::unordered_map<MolKey, vector<size_t>, MolKeyHash> groups;

        for (size_t m = 0; m < components.size(); ++m) {
            auto& atoms_in_mol = components[m];
            size_t n = atoms_in_mol.size();
            std::string name = guess_molecule_name(topo, atoms_in_mol);

            MolKey key{ name, n };
            groups[key].push_back(m);
        }

        Classifier classifier;
        std::vector<internal_segment_desc> comp_desc(components.size());

        for (const auto& kv : groups) {
            const MolKey& key = kv.first;
            const auto& mol_indices = kv.second;
            size_t natom_per_mol = key.natom;
            size_t mol_count = mol_indices.size();

            // One classification result per molblock (group)
            internal_segment_desc group_desc =
                classifier.classify_group(mol_count, natom_per_mol, key.name);

            // Assign same desc to all molecules in this molblock
            for (auto mi : mol_indices) {
                comp_desc[mi] = group_desc;
            }
        }

        // --- Build segments in atom order (this removes contiguity assumption)
        vector<internal_segment_desc> internal_result;
        
        // Start with atom 0
        internal_segment_desc current = comp_desc[comp_id[0]];
        current.size = 1; // number of atoms in this segment

        for (size_t i = 1; i < natoms; ++i) {
            const auto& desc = comp_desc[comp_id[i]];

            // If same type (and optionally same name) -> extend current segment
            // If you only care about type for compression, compare only desc.type.

            if (desc.type == current.type && (desc.group_id == current.group_id || desc.type == internal_segment_desc::segment_type::other)) {
                ++current.size;
            }
            else {
                // Close current segment and start a new one
                internal_result.emplace_back(current);
                current = desc;
                current.size = 1;
            }
        }

        // Push the last segment
        internal_result.emplace_back(current);

        std::size_t mol_id = 1;
        std::size_t wat_id = 1;
        std::size_t other_id = 1;
        std::size_t lip_id = 1;

        vector<segment_desc_t> result;
        result.reserve(internal_result.size());

        for (const auto& x : internal_result)
        {
	        switch (x.type)
	        {
	        case internal_segment_desc::segment_type::other:
                result.emplace_back(mdc::segment_type::other, "other_" + std::to_string(other_id++), x.size);
                break;
            case internal_segment_desc::segment_type::mol:
                result.emplace_back(mdc::segment_type::molecule, "mol_" + std::to_string(mol_id++), x.size);
                break;
            case internal_segment_desc::segment_type::lip:
                result.emplace_back(mdc::segment_type::molecule, "lip_" + std::to_string(lip_id++), x.size);
                break;
            case internal_segment_desc::segment_type::wat:
                result.emplace_back(mdc::segment_type::water, "wat_" + std::to_string(wat_id++), x.size);
                break;
	        default:
                throw std::runtime_error("This should never happen (build desc translate internal segment type)");
	        }
        }
        return result;
    }
    catch (const chemfiles::Error& e) {
        throw std::runtime_error(std::string("Chemfiles error: ") + e.what());
    }
}