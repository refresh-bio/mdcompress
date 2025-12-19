#pragma once

#include "mdc_abstract.h"
#include "engine_decompression.h"

class CMDCDecompress : public CMDCAbstract
{
	CEngineDecompress engine;
	frame_desc_t prev_frame_desc;
	Serializer serializer_desc;

	uint32_t no_frames{};
	std::vector<std::shared_ptr<CEngineDecompress>> segment_engines;
	std::vector<std::vector<uint32_t>> atoms_selected; // atoms_selected[subsegment_id] is a list of selected atoms for subsegment (indexed from 0 per each subsegment)

	std::vector<segment_desc_t> segment_desc;
	std::vector<std::vector<uint32_t>> segment_to_subsegment; // segment_to_subsegment[segement_id] is a list of subsegments ids for segment_id
	std::vector<uint32_t> atom_ids_to_subsegment_ids; // atom_ids_to_subsegment_ids[atom_id] is a subsegment_id

	std::vector<uint32_t> result_atom_id_to_org; // result_atom_id_to_org[result_atom_id] is an original atom_id, this is 0 based!

	packed_frames_t packed_frames;

	std::shared_ptr<CSegmentModel> deserialize_model_molecule(const packed_t& packed);
	std::shared_ptr<CSegmentModel> deserialize_model_water(const packed_t& packed);
	std::shared_ptr<CSegmentModel> deserialize_model_other(const packed_t& packed);

	void allocate_histories_and_engines()
	{
		segment_histories.clear();
		segment_histories.resize(segment_models.size(), CSegmentHistory(max_history_size));

		while (segment_models.size() > segment_engines.size())
			segment_engines.push_back(std::make_shared<CEngineDecompress>());
	}

public:
	CMDCDecompress(uint32_t no_delta_frames, uint32_t max_history_size) : CMDCAbstract(no_delta_frames, max_history_size)
	{}

	bool DeserializeSegmentDesc(const packed_t &packed);

	uint32_t GetNoSubSegments()
	{
		return (uint32_t)subsegment_desc.size();
	}

	const std::vector<segment_desc_t> &GetSegmentDescription() const
	{
		return segment_desc;
	}

	const std::vector<segment_desc_t> &GetSubSegmentDescription() const
	{
		return subsegment_desc;
	}

	bool IsSubsegmentSelected(uint32_t subsegment_id)
	{
		return !atoms_selected[subsegment_id].empty();
	}

	bool DeserializeModel(const uint32_t segment_id, const packed_t& packed);

	bool SelectSegmentsAndAtoms(const std::vector<std::string>& segment_ids, const std::vector<uint32_t>& atom_ids);
	
	bool StartDecompression();

	bool SetPackedFrames(const packed_frames_t& _packed_frames);
	bool DecompressFrame(frame_t& frame);

	bool StartFrame(frame_t& frame, packed_t& packed);
	bool DeserializeAnchorDesc(packed_t& packed);

	const std::vector<uint32_t>& GetAnchorIds() const
	{
		return anchor_ids;
	}

	const std::vector<uint32_t>& GetOriginalAtomIds() const
	{
		return result_atom_id_to_org;
	}

	uint32_t GetNoFrames() const
	{
		return no_frames;
	}
};