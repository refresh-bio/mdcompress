#pragma once

#include <string>
#include "../../libs/libxdrfile/xdrfile_xtc.h"

#include "../common/utils.h"
#include "../common/frame.h"

#if 0

class CXTCReader
{
	XDRFILE* xtc = nullptr;
	int n_atoms{};
	int step{};
	float time{};
	matrix box{};
	float prec{};
	rvec* raw_coords = nullptr;
	bool initialized = false;

	std::vector<segment_desc_t> segment_desc;
	std::vector<coord_t> coords;

public:
	CXTCReader()
	{}

	~CXTCReader()
	{
		if(xtc)
			close();
	}

	bool open(const std::string& fn, const std::vector<segment_desc_t> &_segment_desc)
	{
		if (read_xtc_natoms((char*)fn.c_str(), &n_atoms) != exdrOK)
			return false;

		xtc = xdrfile_open(fn.c_str(), "r");

		if (xtc == nullptr)
			return false;

		raw_coords = new rvec[n_atoms];
		initialized = false;

		segment_desc = _segment_desc;

		uint32_t n_atoms_in_desc = 0;
		for (auto& sd : segment_desc)
			n_atoms_in_desc += sd.size;

		if (n_atoms_in_desc != (uint32_t) n_atoms)
			return false;

		return true;
	}

	bool close()
	{
		if (!xtc)
			return false;

		if (xdrfile_close(xtc) != exdrOK)
			return false;

		xtc = nullptr;

		delete[] raw_coords;
		raw_coords = nullptr;

		return true;
	}

	bool get_frame(frame_t &frame)
	{
		if (xtc == nullptr)
			return false;

		frame.clear();
		frame.desc.n_atoms = n_atoms;

		if (read_xtc(xtc, frame.desc.n_atoms, &frame.desc.step, &frame.desc.time, frame.desc.box, raw_coords, &frame.desc.prec) != exdrOK)
			return false;

		if (!initialized)
		{
			step = frame.desc.step;
			time = frame.desc.time;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box[i][j] = frame.desc.box[i][j];
			prec = frame.desc.prec;

			initialized = true;
		}

		int raw_i = 0;

		for (uint32_t i = 0; i < (uint32_t) segment_desc.size(); ++i)
		{
			uint32_t n_seg_atoms = segment_desc[i].size;

			// Do not load segments marked as NONE
			if (segment_desc[i].type == segment_type_t::none)
			{
				raw_i += n_seg_atoms;
				continue;
			}

			segment_t segment(segment_desc[i].type, segment_desc[i].name, n_seg_atoms, std::vector<coord_t>());

			auto& coords = segment.coords;

			coords.clear();
			coords.reserve(n_seg_atoms);

			for (uint32_t j = 0; j < n_seg_atoms; ++j, ++raw_i)
				coords.emplace_back(iround<int>(raw_coords[raw_i][0] * prec), iround<int>(raw_coords[raw_i][1] * prec), iround<int>(raw_coords[raw_i][2] * prec));

			frame.emplace(std::move(segment));
		}

		return true;
	}
};

#endif
