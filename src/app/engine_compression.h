#pragma once

#include "../common/mdc_abstract.h"
#include "../common/engine.h"

#include "../../libs/refresh/range_coder/lib/rc_utils.h"

class CEngineCompress: public CEngineAbstract
{
	using rc_t = refresh::rc_encoder<refresh::vector_io_stream>;
	CEngineState<rc_t> engine_state;

	using model_prefix_t = CEngineState<rc_t>::model_prefix_t;
	using model_suffix_t = CEngineState<rc_t>::model_suffix_t;

	void rc_encode_int(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix, int64_t val);
	void rc_encode_uint(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix, uint64_t val);

	void compress_other(std::shared_ptr<CSegmentModelOther> model, const segment_t& segment, const CSegmentHistory &history);
	void compress_water(std::shared_ptr<CSegmentModelWater> model, const segment_t& segment, const CSegmentHistory &history);
	void compress_molecule(std::shared_ptr<CSegmentModelMolecule> model, const segment_t& segment, const CSegmentHistory &history);

public:
	CEngineCompress() : CEngineAbstract(),
		engine_state(true)
	{}

	~CEngineCompress()
	{}

	void compress(std::shared_ptr<CSegmentModel> model, const segment_t& segment, const CSegmentHistory &history);
	void open_stream();
	void close_stream(packed_t& packed);
};