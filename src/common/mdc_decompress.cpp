#include <cassert>
#include <set>

#include "mdc_decompress.h"
#include "serializer.h"


// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCDecompress::deserialize_model_molecule(const packed_t& packed)
{
	auto model = std::make_shared<CSegmentModelMolecule>();
	auto& ref_atoms = model->ref_atoms;

	serializer.clear();

	serializer.set_data(packed);
	
	ref_atoms.resize(serializer.load32u());

	for (size_t i = 3; i < ref_atoms.size(); ++i)
	{
		ref_atoms[i].atom_ids[0] = serializer.load8u();
		ref_atoms[i].atom_ids[1] = serializer.load8u();
		ref_atoms[i].atom_ids[2] = serializer.load8u();
		ref_atoms[i].distances[0] = serializer.load32u();
		ref_atoms[i].distances[1] = serializer.load32u();
		ref_atoms[i].distances[2] = serializer.load32u();
		ref_atoms[i].is_left_handed = (bool) serializer.load8u();
	}

	return model;
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCDecompress::deserialize_model_water(const packed_t& packed)
{
	auto model = std::make_shared<CSegmentModelWater>();

	serializer.clear();
	serializer.set_data(packed);

	model->atom_order = serializer.load8u();
	model->dist2_OH = serializer.load64u();
	model->dist2_HH = serializer.load64u();

	return model;
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCDecompress::deserialize_model_other(const packed_t& packed)
{
	auto model = std::make_shared<CSegmentModelOther>();

	return model;
}

// *******************************************************************************
bool CMDCDecompress::DeserializeSegmentDesc(const packed_t &packed)
{
	segment_desc.clear();
	subsegment_desc.clear();
	segment_to_subsegment.clear();	
	serializer.clear();
	atom_ids_to_subsegment_ids.clear();

	serializer.set_data(packed);

	uint32_t n_subsegments = serializer.load32u();

	for (uint32_t i = 0; i < n_subsegments; ++i)
	{
		segment_type_t type = (segment_type_t)serializer.load8u();
		std::string name = serializer.load_str();
		uint32_t n_atoms_in_segment = serializer.load32u();

		subsegment_desc.emplace_back(type, name, n_atoms_in_segment);
	}

	assert(subsegment_desc.size());
	assert(!subsegment_desc.front().name.empty()); //first subsegment must start new segment

	for (uint32_t i = 0 ; i < subsegment_desc.size(); ++i)
	{
		atom_ids_to_subsegment_ids.resize(atom_ids_to_subsegment_ids.size() + subsegment_desc[i].size, i);
		const auto& subseg = subsegment_desc[i];
		//non empty name means begin of new segment
		if (!subseg.name.empty())
		{
			segment_desc.emplace_back(subseg);
			segment_to_subsegment.emplace_back();
			segment_desc.back().size = 0;
		}
		
		segment_desc.back().size += subseg.size;
		segment_to_subsegment.back().push_back(i);
	}

	// Compression features
	compression_features.water_z1 = (compression_features_t::water_z1_t) serializer.load8u();
	compression_features.water_x2 = (compression_features_t::water_x2_t) serializer.load8u();
	compression_features.water_y2 = (compression_features_t::water_y2_t) serializer.load8u();
	compression_features.water_z2 = (compression_features_t::water_z2_t) serializer.load8u();
	compression_features.molecule = (compression_features_t::molecule_t) serializer.load8u();
	compression_features.molecule_orientation = (compression_features_t::molecule_orientation_t) serializer.load8u();
	compression_features.rc_uint32 = (compression_features_t::rc_uint32_t) serializer.load8u();

	return true;
}

// *******************************************************************************
bool CMDCDecompress::DeserializeModel(const uint32_t subsegment_id, const packed_t& packed)
{
	if (subsegment_id >= segment_models.size())
		segment_models.resize(subsegment_id + 1, nullptr);

	switch (subsegment_desc[subsegment_id].type)
	{
	case segment_type_t::molecule:
		segment_models[subsegment_id] = deserialize_model_molecule(packed);
		break;
	case segment_type_t::water:
		segment_models[subsegment_id] = deserialize_model_water(packed);
		break;
	case segment_type_t::other:
		segment_models[subsegment_id] = deserialize_model_other(packed);
		break;
	case segment_type_t::none:
		// Do nothing here
		break;
	default:
		assert(0);
	}

	return true;
}

// *******************************************************************************
bool CMDCDecompress::SelectSegmentsAndAtoms(const std::vector<std::string>& segment_ids, const std::vector<uint32_t>& atom_ids) //atom_ids here is 1 based!
{
	atoms_selected.assign(subsegment_desc.size(), {});

	//if no selection, select all
	if (segment_ids.empty() && atom_ids.empty())
	{
		uint32_t offset = 0;
		uint32_t prev_subseg_id = 0;
		for (uint32_t atom_id = 0 ; atom_id < atom_ids_to_subsegment_ids.size(); ++atom_id)
		{
			const auto subseg_id = atom_ids_to_subsegment_ids[atom_id];
			if (subseg_id != prev_subseg_id)
			{
				offset = atom_id;
				prev_subseg_id = subseg_id;
			}
			atoms_selected[subseg_id].push_back(atom_id - offset);
		}

		result_atom_id_to_org.resize(atom_ids_to_subsegment_ids.size());
		std::iota(result_atom_id_to_org.begin(), result_atom_id_to_org.end(), 1); //start with 1 because in the original notation atom_ids are 1 based
		return true;	
	}

	const auto num_atoms = atom_ids_to_subsegment_ids.size();
	for (auto atom_id : atom_ids)
	{
		assert(atom_id > 0);

		//internally we have zero based!
		--atom_id; //zero based

		if (atom_id >= num_atoms)
			return false;

		const auto subseg_id = atom_ids_to_subsegment_ids[atom_id];

		auto offset = static_cast<uint32_t>(
			std::lower_bound(atom_ids_to_subsegment_ids.begin(), atom_ids_to_subsegment_ids.end(), subseg_id) - atom_ids_to_subsegment_ids.begin());

		atoms_selected[subseg_id].push_back(atom_id - offset);
	}

	std::set<std::string> selected_segments(segment_ids.begin(), segment_ids.end());

	//select all atoms of all subsegments of selected segments
	for (size_t i = 0; i < segment_desc.size(); ++i)
		if (selected_segments.find(segment_desc[i].name) != selected_segments.end())
			for (auto subseg : segment_to_subsegment[i])
			{
				auto first_atom_id = static_cast<uint32_t>(std::lower_bound(atom_ids_to_subsegment_ids.begin(), atom_ids_to_subsegment_ids.end(), subseg) - atom_ids_to_subsegment_ids.begin());
				for (auto atom_id = first_atom_id ; atom_id < atom_ids_to_subsegment_ids.size() && atom_ids_to_subsegment_ids[atom_id] == subseg  ; ++atom_id)
					atoms_selected[subseg].push_back(atom_id - first_atom_id);
			}

	//sort and remove duplicates (very important)
	for (auto& x : atoms_selected)
	{
		std::sort(x.begin(), x.end());
		x.erase(std::unique(x.begin(), x.end()), x.end());
	}

	uint32_t offset = 0;
	
	for (uint32_t ss_id = 0 ; ss_id < subsegment_desc.size() ; ++ss_id)
	{
		for (uint32_t atom_id : atoms_selected[ss_id])
			result_atom_id_to_org.push_back(offset + atom_id + 1); // + 1 because atom_ids in the original notation are 1 based
		offset += subsegment_desc[ss_id].size;
	}
	result_atom_id_to_org.shrink_to_fit();

	return true;
}

// *******************************************************************************
bool CMDCDecompress::StartDecompression()
{
	allocate_histories_and_engines();

	for (auto& se : segment_engines)
		se->set_compression_features(compression_features);

	return true;
}

// *******************************************************************************
bool CMDCDecompress::StartFrame(frame_t& frame, packed_t& packed)
{
	serializer.clear();
	serializer.set_data(packed);

	uint8_t changes = serializer.load8u();

	frame.desc = prev_frame_desc;

	if (changes & 0x01)		frame.desc.n_atoms = serializer.load32i();
	if (changes & 0x02)		frame.desc.step = serializer.load32i();
	if (changes & 0x04)		frame.desc.time = serializer.load32f();
	if (changes & 0x08)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				frame.desc.box[i][j] = serializer.load32f();
	if (changes & 0x10)		frame.desc.prec = serializer.load32f();

	prev_frame_desc = frame.desc;

	return true;
}

// *******************************************************************************
bool CMDCDecompress::SetPackedFrames(const packed_frames_t& _packed_frames)
{
	packed_frames = _packed_frames;

	prev_frame_desc.clear();
	serializer_desc.clear();
	serializer_desc.set_data(packed_frames.desc);

	for (size_t i = 0; i < segment_engines.size(); ++i)
	{
		segment_engines[i]->open_stream(packed_frames.packed[i]);
		segment_histories[i].clear();
	}

	return true;
}

// *******************************************************************************
bool CMDCDecompress::DecompressFrame(frame_t& frame)
{
	frame.segments.resize(subsegment_desc.size());

	// *** Description
	uint8_t changes = serializer_desc.load8u();

	frame.desc = prev_frame_desc;

	if (changes & 0x01)		frame.desc.n_atoms = serializer_desc.load32i();
	if (changes & 0x02)		frame.desc.step = serializer_desc.load32i();
	if (changes & 0x04)		frame.desc.time = serializer_desc.load32f();
	if (changes & 0x08)
	{
		uint16_t non_zeros = serializer_desc.load16u();
		uint16_t m = 1;

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j, m <<= 1)
				if (non_zeros & m)
					frame.desc.box[i][j] = serializer_desc.load64f();
				else
					frame.desc.box[i][j] = 0.0f;
	}

	if (changes & 0x10)		frame.desc.prec = serializer_desc.load32f();

	prev_frame_desc = frame.desc;

	frame.desc.n_atoms = 0;

	// *** Streams
	for (size_t i = 0; i < segment_engines.size(); ++i)
	{
		if (atoms_selected[i].empty())
			continue;

		assert(atoms_selected[i].size() <= subsegment_desc[i].size);

		frame.segments[i].type = subsegment_desc[i].type;
		frame.segments[i].n_atoms = subsegment_desc[i].size; //mkokot_TODO: this value is needed in decompress, but will be adjusted to the real number of selected atoms
		frame.desc.n_atoms += static_cast<int>(atoms_selected[i].size());

//		segment_engines[i]->decompress(segment_models[i], frame.segments[i], segment_histories[i], atoms_selected[i].back());
		segment_engines[i]->decompress(segment_models[i], frame.segments[i], segment_histories[i]);
		segment_histories[i].add(frame.segments[i]);

		if (atoms_selected[i].size() < subsegment_desc[i].size)
		{
			uint32_t write_pos{};
			auto& coords = frame.segments[i].coords;
			for (const auto atom_id : atoms_selected[i]) //atoms_selected[i] must be sorted here
				coords[write_pos++] = coords[atom_id];

			assert(write_pos == atoms_selected[i].size());

			frame.segments[i].coords.resize(atoms_selected[i].size());

			//fix the number of atoms
			frame.segments[i].n_atoms = static_cast<uint32_t>(atoms_selected[i].size());
		}
	}
	return true;
}

// *******************************************************************************
bool CMDCDecompress::DeserializeAnchorDesc(packed_t& packed)
{
	serializer.clear();
	serializer.set_data(packed);

	no_frames = serializer.load32u();
	uint32_t no_anchors = serializer.load32u();

	anchor_ids.clear();

	for (uint32_t i = 0; i < no_anchors; ++i)
		anchor_ids.emplace_back(serializer.load32u());

	return true;
}
