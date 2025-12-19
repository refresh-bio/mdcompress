#pragma once

#include <string>
#include <numeric>
#include <cstring>
#include "../../libs/libxdrfile/xdrfile_xtc.h"

#include "../common/frame.h"


#if 0
class CXTCWriter
{
	XDRFILE* xtc = nullptr;
	int n_atoms{};
	int step{};
	float time{};
	matrix box{};
	float prec{};
	rvec* raw_coords = nullptr;
	bool initialized = false;

public:
	CXTCWriter()
	{}

	~CXTCWriter()
	{
		if (xtc)
			close();
	}

	bool open(const std::string& fn)
	{
		xtc = xdrfile_open(fn.c_str(), "w");

		if (xtc == nullptr)
			return false;

		initialized = false;

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

	//mkokot
	bool add_frames(const mdc::query_result& result)
	{
		if (xtc == nullptr)
			return false;

		for (uint32_t frame_id = 0; frame_id < result.frames.size(); ++frame_id)
		{
			if (!initialized)
			{
				n_atoms = (int)result.frames[frame_id].coords.size();

				raw_coords = new rvec[n_atoms];
				initialized = true;
			}
			else
			{
				if (n_atoms != result.frames[frame_id].coords.size())
				{
					std::cerr << "when writing to xtc the same number of atoms must be used each time\n";
					exit(1);
				}
			}

			for (uint32_t i = 0; i < n_atoms; ++i)
			{
				raw_coords[i][0] = result.frames[frame_id].coords[i].x;
				raw_coords[i][1] = result.frames[frame_id].coords[i].y;
				raw_coords[i][2] = result.frames[frame_id].coords[i].z;
			}


			//all should be the same
			auto& desc = result.frames[frame_id];
			matrix box;
			memcpy(box, desc.box, sizeof(box)); //cannot pass desc.box directly (const)
			if (write_xtc(xtc, n_atoms, desc.step, desc.time, box, raw_coords, desc.prec) != exdrOK)
				return false;
		}
		return true;
	}

	bool add_frame(frame_t& frame)
	{
		if (xtc == nullptr)
			return false;

		if (!initialized)
		{
			n_atoms = 0;
			for (auto& sd : frame.segments)
				n_atoms += (int) sd.coords.size();

			frame.desc.n_atoms = n_atoms;

			raw_coords = new rvec[n_atoms];
			initialized = true;
		}

		uint32_t raw_i = 0;

		for (auto& sd : frame.segments)
		{
			auto& coords = sd.coords;

			for (uint32_t i = 0; i < sd.coords.size(); ++i, ++raw_i)
			{
				raw_coords[raw_i][0] = (float)coords[i].x / frame.desc.prec;
				raw_coords[raw_i][1] = (float)coords[i].y / frame.desc.prec;
				raw_coords[raw_i][2] = (float)coords[i].z / frame.desc.prec;
			}
		}

		if (write_xtc(xtc, n_atoms, frame.desc.step, frame.desc.time, frame.desc.box, raw_coords, frame.desc.prec) != exdrOK)
			return false;

//		frame.clear();

		return true;
	}
};
#endif
