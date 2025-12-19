#pragma once

#include <memory>

#include "mdc_abstract.h"
#include "engine.h"

class CEngineDecompress : public CEngineAbstract
{
	using rc_t = refresh::rc_decoder<refresh::vector_io_stream>;
	CEngineState<rc_t> engine_state;

	using model_prefix_t = CEngineState<rc_t>::model_prefix_t;
	using model_suffix_t = CEngineState<rc_t>::model_suffix_t;

	int64_t rc_decode_uint(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix);
	int64_t rc_decode_int(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix);

	void decompress_other(std::shared_ptr<CSegmentModelOther> model, segment_t& segment, const CSegmentHistory &history);
	void decompress_water(std::shared_ptr<CSegmentModelWater> model, segment_t& segment, const CSegmentHistory &history);
	void decompress_molecule(std::shared_ptr<CSegmentModelMolecule> model, segment_t& segment, const CSegmentHistory &history);

	void decompress_raw(segment_t& segment);

public:
	CEngineDecompress() : CEngineAbstract(),
		engine_state(false)
	{}

	~CEngineDecompress()
	{}

	void decompress(std::shared_ptr<CSegmentModel> model, segment_t& segment, const CSegmentHistory &history);
	void open_stream(const packed_t& packed);
	void close_stream();
};