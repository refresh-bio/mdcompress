#include "mdc_reader.h"
#include <string>
#include <random>
#include <vector>
#include <iostream>
#include <chrono>
#include <fstream>
#include <numeric>

//#define ONLY_XTC_READER
#define USE_CHEMFILES
#include "../app/input_reader.h"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::string_view;
using std::ifstream;

struct Params
{
	string xtc_input_file; //-x
	string mdc_input_file; //-m
	string desc_input_file; //-d
	uint32_t n_frames = 0; //-n, 0 means query all, otherwise query n frames
};


bool parse_params(int argc, char* argv[], Params& params)
{
	if (argc == 1 || (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help")))
	{
		cout << "Usage: " << argv[0] << " -x xtc_input_file -m mdc_input_file -d desc_input_file" << std::endl;
		return false;
	}
	
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]) == "-x")
		{
			if (i + 1 < argc)
			{
				params.xtc_input_file = argv[i + 1];
				i++;
			}
			else
			{
				cerr<< "Error: -x option requires one argument" << std::endl;
				return false;
			}
		}
		else if (string(argv[i]) == "-m")
		{
			if (i + 1 < argc)
			{
				params.mdc_input_file = argv[i + 1];
				i++;
			}
			else
			{
				cerr<< "Error: -m option requires one argument" << std::endl;
				return false;
			}
		}
		else if (string(argv[i]) == "-d")
		{
			if (i + 1 < argc)
			{
				params.desc_input_file = argv[i + 1];
				i++;
			}
			else
			{
				cerr<< "Error: -d option requires one argument" << std::endl;
				return false;
			}
		}
		else if (string(argv[i]) == "-n")
		{
			if (i + 1 < argc)
			{
				params.n_frames = std::stoi(argv[i + 1]);
				i++;
			}
			else
			{
				cerr<< "Error: -n option requires one argument" << std::endl;
				return false;
			}
		}
		else
		{
			cerr<< "Error: unknown option " << argv[i] << std::endl;
			return false;
		}
	}

	if (params.xtc_input_file.empty())
	{
		cerr<< "Error: -x is required" << std::endl;
		return false;
	}

	if (params.mdc_input_file.empty())
	{
		cerr<< "Error: -m is required" << std::endl;
		return false;
	}

	if (params.desc_input_file.empty())
	{
		cerr<< "Error: -d is required" << std::endl;
		return false;
	}

	return true;
}

// *******************************************************************************
using segment_type_t = mdc::segment_type;
using segment_desc_t = mdc::segment_desc_t;
bool load_segment_description_file(const std::string& path, vector<segment_desc_t>& segment_desc)
{
	ifstream ifs(path);

	if (!ifs)
	{
		cerr << "Cannot open file: " << path  << endl;
		return false;
	}

	segment_desc.clear();

	string s_type;
	segment_type_t type;
	string name;
	uint32_t n_atoms;

	while (ifs >> s_type >> name >> n_atoms)
	{
		if (s_type == std::string("MOLECULE") || s_type == std::string("MOL"))
			type = segment_type_t::molecule;
		else if (s_type == std::string("WATER") || s_type == std::string("WAT"))
			type = segment_type_t::water;
		else if (s_type == std::string("OTHER") || s_type == std::string("OTH"))
			type = segment_type_t::other;
		else if (s_type == std::string("NONE") || s_type == std::string("NON"))
			type = segment_type_t::none;
		else
		{
			cerr << "Unknown segment type: " << s_type << endl;
			return false;
		}

		segment_desc.emplace_back(type, name, n_atoms);
	}

	return true;
}

bool load_xtc_file(const std::string& path, const vector<segment_desc_t>& segment_desc, vector<frame_t>& frames)
{
#ifdef USE_CHEMFILES
	CTrajReader traj_reader(2);
	if (!traj_reader.open(path, segment_desc))
	{
		std::cerr << "Cannot open: " << path << std::endl;
		return false;
	}

	frame_t frame;

	while (traj_reader.get_frame(frame))
		frames.push_back(std::move(frame));

	traj_reader.close();

	return true;
#else
	CXTCReader xtc_reader;
	if (!xtc_reader.open(path, segment_desc))
	{
		cerr << "Cannot open: " << path << endl;
		return false;
	}

	frame_t frame;
	frames.clear();
	for (uint32_t frame_id = 0; xtc_reader.get_frame(frame); ++frame_id)
	{
		frames.emplace_back(std::move(frame));
	}
	return true;
#endif
}

bool operator==(const segment_desc_t& lhs, const segment_desc_t& rhs)
{
	if (lhs.size != rhs.size)
		return false;
	if (lhs.name != rhs.name)
		return false;
	if (lhs.type != rhs.type)
		return false;
	return true;
}

bool compare_frames(const frame_t& lhs, const mdc::frame& rhs, const std::vector<uint32_t>& to_org_atom_ids)
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (lhs.desc.box[i][j] != rhs.box[i][j])
				return false;

//	if (memcmp(lhs.desc.box, rhs.box, sizeof(lhs.desc.box)) != 0)
//		return false;

	//may be different if we ask for specific atoms
	//if (lhs.desc.n_atoms != rhs.coords.size())
	//	return false;

	if (lhs.desc.time != rhs.time)
		return false;

	if(lhs.desc.prec != rhs.prec)
		return false;

	if(lhs.desc.step != rhs.step)
		return false;

	for (uint32_t res_atom_ids = 0; res_atom_ids < rhs.coords.size(); ++res_atom_ids)
	{
		uint32_t org_atom_id = to_org_atom_ids[res_atom_ids] - 1;

		//I'm just looping but this is to be rewritten to have lookup table for atom ids to segments
		uint32_t cur_atom_id = 0;
		bool B = false;
		for (const auto& seg : lhs.segments)
		{
			for (const auto& atom : seg.coords)
			{
				if (cur_atom_id == org_atom_id)
				{
					if (atom.x / rhs.prec != rhs.coords[res_atom_ids].x)
						return false;
					if (atom.y / rhs.prec != rhs.coords[res_atom_ids].y)
						return false;
					if (atom.z / rhs.prec != rhs.coords[res_atom_ids].z)
						return false;
					B = true;
					break;
				}
				cur_atom_id++;
			}
			if (B)
				break;
		}
	}
	return true;
}

bool check(
	const vector<frame_t>& frames,
	const vector<segment_desc_t>& segment_desc,
	mdc::reader& reader,
	uint32_t query_for_n_frames, //0 - query for all
	const std::vector<std::string>& segment_ids,
	const std::vector<uint32_t>& atom_ids)
{
	auto query_engine_ptr = reader.get_query_engine(segment_ids, atom_ids);
	auto& query_engine = *query_engine_ptr;
	if (const auto& err = query_engine.get_current_error(); !err.empty())
	{
		std::cerr << err;
		exit(1);
	}
	uint32_t no_segments = reader.get_segments().size();
	uint32_t no_frames = reader.get_no_frames();

	if (segment_desc.size() != no_segments)
	{
		cerr << "Error: segment description file does not match the number of segments in the MDC file" << endl;
		exit(1); 
	}

	if (no_frames != frames.size())
	{
		cerr << "Error: number of frames in the XTC file does not match the number of frames in the MDC file" << endl;
		exit(1);
	}

	//auto a = segment_desc.front() == query_engine.get_segments().front();
	if(segment_desc.size() != reader.get_segments().size())
	{
		cerr << "Error: segment description file does not match the number of segments in the MDC file" << endl;
		exit(1);
	}
	for(uint32_t i = 0; i < segment_desc.size(); ++i)
	{
		if (!(segment_desc[i] == reader.get_segments()[i]))
		{
			cerr << "Error: segment description file does not match segments in the MDC file" << endl;
			exit(1);
		}
	}

	vector queried_atom_ids(atom_ids);
	uint32_t q_id = 0;
	for (uint32_t i = 0; i < segment_desc.size(); ++i)
	{
		//if querying for fiven segment desc
		if (std::find(segment_ids.begin(), segment_ids.end(), segment_desc[i].name) != segment_ids.end())
			for (uint32_t j = 0; j < segment_desc[i].size; ++j)
				queried_atom_ids.push_back(q_id++);
		else
			q_id += segment_desc[i].size;
	}
	std::sort(queried_atom_ids.begin(), queried_atom_ids.end());
	queried_atom_ids.erase(std::unique(queried_atom_ids.begin(), queried_atom_ids.end()), queried_atom_ids.end());

	if (queried_atom_ids == query_engine.get_original_atom_ids())
	{
		std::cerr << "Error: queried atom ids do not match the original atom ids" << std::endl;
		return false;
	}
	
	mdc::query_result result;

	std::vector<uint32_t> frame_ids(no_frames);
	std::iota(frame_ids.begin(), frame_ids.end(), 0ul);
	if (query_for_n_frames > 0 && query_for_n_frames < no_frames)
		frame_ids.resize(query_for_n_frames);
		
	std::uniform_int_distribution<uint32_t> dist_no_frames(1, 10); //query for a a random number of frames at once
	std::shuffle(frame_ids.begin(), frame_ids.end(), std::mt19937{});
	std::mt19937 mt;
	for (uint32_t frame_id_id = 0; frame_id_id < no_frames; )
	{
		auto n_frames = dist_no_frames(mt);
		auto end = frame_id_id + n_frames;
		if (end >  no_frames)
			end = no_frames;

		std::span frame_ids_span(frame_ids.data() + frame_id_id, end - frame_id_id);

		frame_id_id = end;
		
		if (!query_engine.query(frame_ids_span, result))
		{
			std::cerr << "Error: " << query_engine.get_current_error();
			exit(1);
		}

		if (result.frames.size() != frame_ids_span.size())
		{
			std::cerr << "Error: query result does not contain the expected number of frames" << std::endl;
			exit(1);
		}

		uint32_t res_i{};
		for (uint32_t frame_id : frame_ids_span)
		{
			auto& result_frame = result.frames[res_i++];

			auto& frame_pattern = frames[frame_id];

			if (!compare_frames(frame_pattern, result_frame, query_engine.get_original_atom_ids()))
			{
				std::cerr << "Error: frame " << frame_id << " does not match" << std::endl;
				return false;
			}
		}
	}
	return true;
}

int main(int argc, char* argv[])
{
	Params params;
	if (!parse_params(argc, argv, params))
		return 1;

	cout << "xtc_input_file: " << params.xtc_input_file << std::endl;
	cout << "mdc_input_file: " << params.mdc_input_file << std::endl;
	cout << "desc_input_file: " << params.desc_input_file << std::endl;
	cout << "n_frames: " << params.n_frames << std::endl;

	vector<segment_desc_t> segment_desc;
	if (!load_segment_description_file(params.desc_input_file, segment_desc))
		return 1;

	vector<frame_t> frames;

	std::cerr << "Reading xtc file...";
	auto start_time = std::chrono::high_resolution_clock::now();
	if (!load_xtc_file(params.xtc_input_file, segment_desc, frames))
		return 1;
	auto end_time = std::chrono::high_resolution_clock::now();

	std::cerr << "\nDone. Time: " << std::chrono::duration<double>(end_time - start_time).count() << "\n";

	mdc::reader mdc_reader(params.mdc_input_file);

	mdc::reader reader(params.mdc_input_file);
	if (const auto& err = reader.get_current_error(); !err.empty())
	{
		std::cerr << err;
		return 1;
	}

	
	vector<vector<string>> segment_ids_cfg;
	vector<vector<uint32_t>> atom_ids_cfg;
	
	uint32_t n_atoms{};
	for (const auto& x : segment_desc)
		n_atoms += x.size;

	cout << "n_atoms: " << n_atoms << std::endl;
	
	//query for all frames and atoms
	segment_ids_cfg.emplace_back(); 
	atom_ids_cfg.emplace_back();

	std::vector<uint32_t> all_atom_ids(n_atoms);
	std::iota(all_atom_ids.begin(), all_atom_ids.end(), 1ull); //indexed from 1

	atom_ids_cfg.emplace_back(all_atom_ids); //all atom ids in order
	
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + n_atoms / 2); //half
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + n_atoms / 4); // quarter
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + 10); // just 10 atoms
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + 1); // a single

	//the same but for shuffled
	std::shuffle(all_atom_ids.begin(), all_atom_ids.end(), std::mt19937{});
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + n_atoms / 2); //half
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + n_atoms / 4); // quarter
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + 10); // just 10 atoms
	atom_ids_cfg.emplace_back(all_atom_ids.begin(), all_atom_ids.begin() + 1); // a single

	vector<string> all_segment_names;
	all_segment_names.reserve(segment_desc.size());
	for (const auto& x : segment_desc)
		all_segment_names.push_back(x.name);

	segment_ids_cfg.emplace_back(all_segment_names); //all segments

	for (const auto& x : segment_desc)
		segment_ids_cfg.emplace_back(std::vector<std::string>{x.name}); //single segment

	if (all_segment_names.size() == 2)
	{
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[0], all_segment_names[1]});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[1], all_segment_names[0]});
	}

	if (all_segment_names.size() == 3)
	{
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[0], all_segment_names[2]});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[2], all_segment_names[0]});

		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[0], all_segment_names[1], all_segment_names[2]});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[0], all_segment_names[2], all_segment_names[1]});

		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[1], all_segment_names[0], all_segment_names[2]});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[1], all_segment_names[2], all_segment_names[0]});

		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[2], all_segment_names[0], all_segment_names[1]});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[2], all_segment_names[1], all_segment_names[0]});

	}
	else if (all_segment_names.size() > 3)
	{
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names[0], all_segment_names.back()});
		segment_ids_cfg.emplace_back(std::vector<std::string>{all_segment_names.back(), all_segment_names.front()});
	}

	//sort segment ids cfg and atoms ids cfg such that we start with smallest vectors
	auto cmp = [](const auto& lhs, const auto& rhs)
	{
		//treat empty as largest possible, becuse empty vec means take all
		if (lhs.size() == 0 && rhs.size() == 0)
			return false;
		if (lhs.size() == 0)
			return false;
		if (rhs.size() == 0)
			return true;

		return lhs.size() < rhs.size();
	};

	std::sort(segment_ids_cfg.begin(), segment_ids_cfg.end(), cmp);
	std::sort(atom_ids_cfg.begin(), atom_ids_cfg.end(), cmp);
	
	for (const auto& segment_ids : segment_ids_cfg)
	{
		std::cerr << "query for segments: ";
		for (auto& x : segment_ids)
			std::cerr << x << " ";
		std::cerr << "\n";

		for (const auto& atom_ids : atom_ids_cfg)
		{
			std::cerr << "\t#atoms: " << atom_ids.size();
			if (atom_ids.size() < 30)
			{
				std::cerr << ", atom ids: ";
				for(auto x : atom_ids)
					std::cerr << x << " ";
			}
			std::cerr << "\n";

			if (!check(frames, segment_desc, reader, params.n_frames, segment_ids, atom_ids))
			{
				//can print more info here
				std::cerr << "Error: check failed" << std::endl;
				return 1;
			}
		}
	}
		

	return 0;
}