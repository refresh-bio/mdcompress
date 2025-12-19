#include <algorithm>
#include <memory>

#include "mdc_compress.h"

// *******************************************************************************
void CMDCCompressWorker::serialize_frame_desc(const frame_desc_t& desc)
{
	uint8_t changes = 0;

	if (prev_frame_desc.n_atoms != desc.n_atoms)
	{
		changes |= 0x01;
		prev_frame_desc.n_atoms = desc.n_atoms;
	}

	if (prev_frame_desc.step != desc.step)
	{
		changes |= 0x02;
		prev_frame_desc.step = desc.step;
	}

	if (prev_frame_desc.time != desc.time)
	{
		changes |= 0x04;
		prev_frame_desc.time = desc.time;
	}

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (prev_frame_desc.box[i][j] != desc.box[i][j])
			{
				changes |= 0x08;
				prev_frame_desc.box[i][j] = desc.box[i][j];
			}

	if (prev_frame_desc.prec != desc.prec)
	{
		changes |= 0x10;
		prev_frame_desc.prec = desc.prec;
	}

	serializer_desc.append8u(changes);

	if (changes & 0x01)		serializer_desc.append32i(desc.n_atoms);
	if (changes & 0x02)		serializer_desc.append32i(desc.step);
	if (changes & 0x04)		serializer_desc.append32f(desc.time);
	if (changes & 0x08)
	{
		uint16_t non_zeros = 0;
		uint16_t m = 1;

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j, m <<= 1)
				if (desc.box[i][j] != 0.0f)
					non_zeros |= m;

		serializer_desc.append16u(non_zeros);

		m = 1;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j, m <<= 1)
				if (non_zeros & m)
					serializer_desc.append64f(desc.box[i][j]);
	}

	if (changes & 0x10)		serializer_desc.append32f(desc.prec);
}

// *******************************************************************************
void CMDCCompressWorker::StartCompression(
	const std::vector<std::shared_ptr<CSegmentModel>>& segment_models,
	uint32_t max_history_size,
	const compression_features_t& compression_features)
{
	this->segment_models = &segment_models;
	allocate_histories_and_engines(segment_models, max_history_size);
	for (auto& se : segment_engines)
		se->set_compression_features(compression_features);
}

// *******************************************************************************
void CMDCCompressWorker::HandleSegment(
	uint32_t segment_id,
	const std::vector<frame_t>& batch,
	packed_t& packed)
{
	const std::vector<std::shared_ptr<CSegmentModel>>& segment_models = *this->segment_models;
	segment_engines[segment_id]->open_stream();
	segment_histories[segment_id].clear();

	const auto no_frames_in_batch = static_cast<uint32_t>(batch.size());
	for (uint32_t frame_id = 0; frame_id < no_frames_in_batch; ++frame_id)
	{
		const auto& cur_frame = batch[frame_id];
		segment_engines[segment_id]->compress(segment_models[segment_id], cur_frame.segments[segment_id], segment_histories[segment_id]);
		segment_histories[segment_id].add(cur_frame.segments[segment_id]);
	}

	segment_engines[segment_id]->close_stream(packed);
}

// *******************************************************************************
void CMDCCompressWorker::HandleDesc(const std::vector<frame_t>& batch,
	packed_t& packed)
{
	prev_frame_desc.clear();
	serializer_desc.clear();
	for (const auto& frame : batch)
		serialize_frame_desc(frame.desc);

	packed = serializer_desc.get_data();
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCCompress::build_model_molecule(const segment_t& segment, const uint32_t max_dist_in_model)
{
	auto model = std::make_shared<CSegmentModelMolecule>();

	auto& ref_atoms = model->ref_atoms;

	ref_atoms.reserve(segment.n_atoms);

	std::vector<std::pair<uint64_t, uint64_t>> d2;

	d2.reserve(max_dist_in_model);

	ref_atoms.resize(std::min<uint32_t>(3, segment.n_atoms));

	for (uint32_t i = 3; i < segment.n_atoms; ++i)
	{
		d2.clear();

		for (uint32_t j = 1; j <= max_dist_in_model && j <= i; ++j)
			d2.emplace_back(dist2_3d(segment.coords[i], segment.coords[i - j]), j);

		partial_sort(d2.begin(), d2.begin() + 3, d2.end());
		ref_atoms.emplace_back((uint32_t) d2[0].second, (uint32_t) d2[1].second, (uint32_t) d2[2].second, 
			uround<uint32_t>(sqrt(d2[0].first)), uround<uint32_t>(sqrt(d2[1].first)), uround<uint32_t>(sqrt(d2[2].first)));
	}

	return model;
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCCompress::build_model_molecule(const std::vector<segment_t>& segments, const uint32_t max_dist_in_model)
{
	auto model = std::make_shared<CSegmentModelMolecule>();

	auto& ref_atoms = model->ref_atoms;

	ref_atoms.reserve(segments.front().n_atoms);

	std::vector<std::pair<uint64_t, uint64_t>> d2;

	d2.reserve(max_dist_in_model);

	ref_atoms.resize(std::min<uint32_t>(3, segments.front().n_atoms));

	assert(max_dist_in_model >= 3);

	for (uint32_t i = 3; i < segments.front().n_atoms; ++i)
	{
		d2.clear();
		for (uint32_t j = 1; j <= max_dist_in_model && j <= i; ++j)
			d2.emplace_back(0, j);

		for(uint32_t k = 0; k < (uint32_t) segments.size(); ++k)
			for (uint32_t j = 1; j <= max_dist_in_model && j <= i; ++j)
				d2[j-1].first += dist_3d(segments[k].coords[i], segments[k].coords[i - j]);

		partial_sort(d2.begin(), d2.begin() + 3, d2.end());
		ref_atoms.emplace_back((uint32_t) d2[0].second, (uint32_t) d2[1].second, (uint32_t) d2[2].second, 
			uround_div<uint32_t>(d2[0].first, uint64_t(segments.size())), uround_div<uint32_t>(d2[1].first, uint64_t(segments.size())), uround_div<uint32_t>(d2[2].first, uint64_t(segments.size())));
	}

	return model;
}

// *******************************************************************************
void CMDCCompress::add_orientation(const uint32_t subsegment_id, const segment_t& segment)
{
	auto& ref_atoms = dynamic_pointer_cast<CSegmentModelMolecule>(segment_models[subsegment_id])->ref_atoms;

	for (uint32_t i = 3; i < segment.n_atoms; ++i)
	{
		auto& cri = ref_atoms[i].atom_ids;

		ref_atoms[i].is_left_handed = engine.is_left_handed(segment.coords[i - cri[0]], segment.coords[i - cri[1]], segment.coords[i - cri[2]], segment.coords[i]);
	}
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCCompress::build_model_water(const segment_t& segment)
{
	auto model = std::make_shared<CSegmentModelWater>();

	uint64_t dist2_01 = 0;
	uint64_t dist2_02 = 0;
	uint64_t dist2_12 = 0;

	std::vector<uint64_t> v_dist2_01;
	std::vector<uint64_t> v_dist2_02;
	std::vector<uint64_t> v_dist2_12;

	auto& coords = segment.coords;

	for (uint32_t i = 0; i < coords.size(); i += 3)
	{
		v_dist2_01.emplace_back(dist2_3d(coords[i], coords[i + 1]));
		v_dist2_02.emplace_back(dist2_3d(coords[i], coords[i + 2]));
		v_dist2_12.emplace_back(dist2_3d(coords[i + 1], coords[i + 2]));
	}

	nth_element(v_dist2_01.begin(), v_dist2_01.begin() + v_dist2_01.size() / 2, v_dist2_01.end());
	nth_element(v_dist2_02.begin(), v_dist2_02.begin() + v_dist2_02.size() / 2, v_dist2_02.end());
	nth_element(v_dist2_12.begin(), v_dist2_12.begin() + v_dist2_12.size() / 2, v_dist2_12.end());

	dist2_01 = *(v_dist2_01.begin() + v_dist2_01.size() / 2);
	dist2_02 = *(v_dist2_02.begin() + v_dist2_02.size() / 2);
	dist2_12 = *(v_dist2_12.begin() + v_dist2_12.size() / 2);

	if (dist2_01 > dist2_02 && dist2_01 > dist2_12)
	{
		model->atom_order = 2;		// HHO
		model->dist2_HH = dist2_01;
		model->dist2_OH = (dist2_02 + dist2_12) / 2;
	}
	else if (dist2_02 > dist2_01 && dist2_02 > dist2_12)
	{
		model->atom_order = 1;		// HOH
		model->dist2_HH = dist2_02;
		model->dist2_OH = (dist2_01 + dist2_12) / 2;
	}
	else 
	{
		model->atom_order = 0;		// OHH
		model->dist2_HH = dist2_12;
		model->dist2_OH = (dist2_01 + dist2_02) / 2;
	}

	return model;
}

// *******************************************************************************
std::shared_ptr<CSegmentModel> CMDCCompress::build_model_other(const segment_t& segment)
{
	return std::make_shared<CSegmentModelOther>();
}

// *******************************************************************************
packed_t CMDCCompress::serialize_model_molecule(std::shared_ptr<CSegmentModelMolecule> model)
{
	serializer.clear();

	auto& ref_atoms = model->ref_atoms;

	serializer.append32u((uint32_t) ref_atoms.size());

	for (size_t i = 3; i < ref_atoms.size(); ++i)
	{
		serializer.append8u(ref_atoms[i].atom_ids[0]);
		serializer.append8u(ref_atoms[i].atom_ids[1]);
		serializer.append8u(ref_atoms[i].atom_ids[2]);
		serializer.append32u(ref_atoms[i].distances[0]);
		serializer.append32u(ref_atoms[i].distances[1]);
		serializer.append32u(ref_atoms[i].distances[2]);
		serializer.append8u(uint8_t(ref_atoms[i].is_left_handed));
	}

	return serializer.get_data();
}

// *******************************************************************************
packed_t CMDCCompress::serialize_model_water(std::shared_ptr<CSegmentModelWater> model)
{
	serializer.clear();

	serializer.append8u(model->atom_order);
	serializer.append64u(model->dist2_OH);
	serializer.append64u(model->dist2_HH);

	return serializer.get_data();
}

// *******************************************************************************
packed_t CMDCCompress::serialize_model_other(std::shared_ptr<CSegmentModelOther> model)
{
	return packed_t();
}

// *******************************************************************************
void CMDCCompress::SerializeSegmentDesc(packed_t& packed)
{
	serializer.clear();

	serializer.append32u((uint32_t)subsegment_desc.size());

	for (auto& sd : subsegment_desc)
	{
		serializer.append8u((uint8_t)sd.type);
		serializer.append_str(sd.name);
		serializer.append32u(sd.size);
	}

	// Compression features
	serializer.append8u((uint8_t)compression_features.water_z1);
	serializer.append8u((uint8_t)compression_features.water_x2);
	serializer.append8u((uint8_t)compression_features.water_y2);
	serializer.append8u((uint8_t)compression_features.water_z2);
	serializer.append8u((uint8_t)compression_features.molecule);
	serializer.append8u((uint8_t)compression_features.molecule_orientation);
	serializer.append8u((uint8_t)compression_features.rc_uint32);

	packed = serializer.get_data();
}

// *******************************************************************************
void CMDCCompress::BuildModel(const uint32_t subsegment_id, const segment_t& segment, const CParams &params)
{
	if (subsegment_id >= segment_models.size())
	{
		segment_models.resize(subsegment_id + 1, nullptr);
		subsegment_desc.resize(subsegment_id + 1);
	}

	subsegment_desc[subsegment_id].type = segment.type;
	subsegment_desc[subsegment_id].name = segment.name;
	subsegment_desc[subsegment_id].size = segment.n_atoms;

	switch (segment.type)
	{
	case segment_type_t::molecule:
		segment_models[subsegment_id] = build_model_molecule(segment, params.max_dist_in_model());
		break;
	case segment_type_t::water:
		segment_models[subsegment_id] = build_model_water(segment);
		break;
	case segment_type_t::other:
		segment_models[subsegment_id] = build_model_other(segment);
		break;
	case segment_type_t::none:
		// Do nothing here
		break;
	default:
		assert(0);
	}
}

// *******************************************************************************
void CMDCCompress::BuildModel(const uint32_t subsegment_id, const std::vector<segment_t>& segments, const CParams& params)
{
	if (subsegment_id >= segment_models.size())
	{
		segment_models.resize(subsegment_id + 1, nullptr);
		subsegment_desc.resize(subsegment_id + 1);
	}

	subsegment_desc[subsegment_id].type = segments.front().type;
	subsegment_desc[subsegment_id].name = segments.front().name;
	subsegment_desc[subsegment_id].size = segments.front().n_atoms;

	switch (segments.front().type)
	{
	case segment_type_t::molecule:
		segment_models[subsegment_id] = build_model_molecule(segments, params.max_dist_in_model());
		break;
	case segment_type_t::water:
		segment_models[subsegment_id] = build_model_water(segments.front());
		break;
	case segment_type_t::other:
		segment_models[subsegment_id] = build_model_other(segments.front());
		break;
	case segment_type_t::none:
		// Do nothing here
		break;
	default:
		assert(0);
	}
}

// *******************************************************************************
void CMDCCompress::AddOrientation(const uint32_t subsegment_id, const segment_t& segment, const CParams& params)
{
	switch (segment.type)
	{
	case segment_type_t::molecule:
		add_orientation(subsegment_id, segment);
		break;
	case segment_type_t::water: [[fallthrough]];
	case segment_type_t::other: [[fallthrough]];
	case segment_type_t::none:
		break;
	default:
		assert(0);
	}
}

// *******************************************************************************
bool CMDCCompress::SerializeModel(const uint32_t segment_id, packed_t& packed)
{
	switch (subsegment_desc[segment_id].type)
	{
	case segment_type_t::molecule:
		packed = serialize_model_molecule(dynamic_pointer_cast<CSegmentModelMolecule>(segment_models[segment_id]));
		break;
	case segment_type_t::water:
		packed = serialize_model_water(dynamic_pointer_cast<CSegmentModelWater>(segment_models[segment_id]));
		break;
	case segment_type_t::other:
		packed = serialize_model_other(dynamic_pointer_cast<CSegmentModelOther>(segment_models[segment_id]));
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
bool CMDCCompress::StartCompression()
{
	no_frames = 0;
	priority = 0;

	return true;
}

// *******************************************************************************
void CMDCCompress::SplitBatch(std::vector<frame_t>& batch, std::vector<compress_task_desc_t>& generated_tasks)
{
	assert(!batch.empty());

	generated_tasks.clear();

	anchor_ids.push_back(no_frames);
	no_frames += batch.size();

	std::shared_ptr<std::vector<frame_t>> batch_sp = std::make_shared<std::vector<frame_t>>(std::move(batch));
	for (size_t segment_id = 0; segment_id < batch_sp->front().segments.size(); ++segment_id)
	{
		generated_tasks.emplace_back();
		compress_task_desc_t& task = generated_tasks.back();

		task.priority = priority++;
		task.batch = batch_sp;
		task.work_type = compress_task_desc_t::work_type_t::segment;
		task.segment_id = segment_id;
	}

	generated_tasks.emplace_back();
	compress_task_desc_t& task = generated_tasks.back();
	task.priority = priority++;
	task.batch = batch_sp;
	task.work_type = compress_task_desc_t::work_type_t::desc;
	task.segment_id = -1;
}
// *******************************************************************************
bool CMDCCompress::Finalize(packed_t& anchors_desc)
{
	Serializer serializer;

	serializer.append32u(no_frames);
	serializer.append32u((uint32_t)anchor_ids.size());
	for (auto x : anchor_ids)
		serializer.append32u(x);

	anchors_desc = serializer.get_data();

	return true;
}
