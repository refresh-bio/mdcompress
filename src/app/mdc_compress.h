#pragma once

#include <mimalloc/atomic.h>

#include "params.h"
#include "../common/mdc_abstract.h"
#include "engine_compression.h"

class CMDCCompressWorker
{
	std::vector<std::shared_ptr<CEngineCompress>> segment_engines;
	std::vector<CSegmentHistory> segment_histories;

	Serializer serializer_desc;
	frame_desc_t prev_frame_desc;

	void allocate_histories_and_engines(const std::vector<std::shared_ptr<CSegmentModel>>& segment_models, uint32_t max_history_size)
	{
		segment_histories.clear();
		segment_histories.resize(segment_models.size(), CSegmentHistory(max_history_size));

		while (segment_models.size() > segment_engines.size())
			segment_engines.push_back(std::make_shared<CEngineCompress>());
	}
	void serialize_frame_desc(const frame_desc_t& desc);
	const std::vector<std::shared_ptr<CSegmentModel>>* segment_models{};
public:
	void StartCompression(
		const std::vector<std::shared_ptr<CSegmentModel>>& segment_models,
		uint32_t max_history_size,
		const compression_features_t& compression_features);

	void HandleSegment(
		uint32_t segment_id,
		const std::vector<frame_t>& batch,
		packed_t& packed);

	void HandleDesc(const std::vector<frame_t>& batch,
		packed_t& packed);
};


struct compress_task_desc_t
{
	size_t priority;
	std::shared_ptr<std::vector<frame_t>> batch;

	enum class work_type_t { segment, desc } work_type;
	size_t segment_id;
};



class CMDCCompress : public CMDCAbstract
{
	CEngineCompress engine;

	uint32_t no_frames{};

	std::shared_ptr<CSegmentModel> build_model_molecule(const segment_t& segment, const uint32_t max_dist_in_model);
	std::shared_ptr<CSegmentModel> build_model_molecule(const std::vector<segment_t>& segments, const uint32_t max_dist_in_model);
	std::shared_ptr<CSegmentModel> build_model_water(const segment_t& segment);
	std::shared_ptr<CSegmentModel> build_model_other(const segment_t& segment);

	void add_orientation(const uint32_t subsegment_id, const segment_t& segment);

	packed_t serialize_model_molecule(std::shared_ptr<CSegmentModelMolecule> model);
	packed_t serialize_model_water(std::shared_ptr<CSegmentModelWater> model);
	packed_t serialize_model_other(std::shared_ptr<CSegmentModelOther> model);

	size_t priority{};
public:
	CMDCCompress(uint32_t no_delta_frames, uint32_t max_history_size, compression_features_t compression_features) :
		CMDCAbstract(no_delta_frames, max_history_size, compression_features)
	{}

	~CMDCCompress()
	{}

	void InitWorker(CMDCCompressWorker& worker)
	{
		worker.StartCompression(segment_models, max_history_size, compression_features);
	}

	void SerializeSegmentDesc(packed_t& packed);

	void BuildModel(const uint32_t subsegment_id, const segment_t& segment, const CParams &params);
	void BuildModel(const uint32_t subsegment_id, const std::vector<segment_t>& segment, const CParams &params);
	void AddOrientation(const uint32_t subsegment_id, const segment_t& segment, const CParams& params);
	bool SerializeModel(const uint32_t segment_id, packed_t& packed);

	bool StartCompression();
	bool Finalize(packed_t &anchors_desc);

	void SplitBatch(std::vector<frame_t>& batch, std::vector<compress_task_desc_t>& generated_tasks);

	size_t GetNoFrames() const
	{
		return no_frames;
	}

	size_t GetNoAtomsInSubSegment(size_t i) const
	{
		return subsegment_desc[i].size;
	}
};
