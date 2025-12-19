#pragma once

#include <cassert>
#include "params.h"
#include "xtc_reader.h"
#include "xtc_writer.h"
#include "mdc_compress.h"
#include "../common/mdc_decompress.h"

#include "../libs/refresh/archive/lib/archive_input.h"
#include "../libs/refresh/archive/lib/archive_output.h"
#include "../libs/refresh/parallel_queues/lib/parallel-queues.h"

class CApplication
{
	CParams params;

	std::vector<segment_desc_t> segment_desc;

	bool load_segment_description_file();

	bool load_desc_from_trajectory_file();

	void reduce_segment_desc_to_molecules_only();

	void expand_segments();

	bool compress();
	bool decompress() const;
	bool select() const;
	bool convert()const;
	bool make_desc()const;
	std::vector<uint32_t> make_frame_ids(uint32_t tot_frames) const;
	bool info();
	bool to_fp32();

	std::pair<uint32_t, uint32_t> locate_frame(const std::vector<uint32_t>& anchor_ids, const uint32_t frame_id) const
	{
		auto p = upper_bound(anchor_ids.begin(), anchor_ids.end(), frame_id) - 1;
		return std::make_pair(p - anchor_ids.begin(), frame_id - *p);
	}

	const char* segment_type_t_to_str(const segment_type_t type) const;

public:
	CApplication()
	{}

	bool run(const CParams& _params)
	{
		params = _params;

		switch (params.mode)
		{
		case CParams::mode_t::compress:
			return compress();
		case CParams::mode_t::decompress:
			return decompress();
		case CParams::mode_t::select:
			return select();
		case CParams::mode_t::make_desc:
			return make_desc();
		case CParams::mode_t::convert:
			return convert();
		case CParams::mode_t::info:
			return info();
		case CParams::mode_t::to_fp32:
			return to_fp32();
		default:
			assert(0);			// Never should be here
		}

		return false;
	}
};
