#pragma once

#include <string>
#include <chemfiles.hpp>

#include "../common/utils.h"
#include "../common/frame.h"

#if 0

class CTrajReader
{
	chemfiles::Trajectory* traj = nullptr;
	chemfiles::Frame ch_frame;

	bool ch_frame_buffered;
	int n_atoms{};
	int step{};
	float time{};
	matrix box{};
	float prec{};
	bool initialized = false;

	std::vector<segment_desc_t> segment_desc;
	std::vector<coord_t> coords;

public:
	CTrajReader()
	{}

	~CTrajReader()
	{
		if(traj)
			close();
	}

	bool open(const std::string& fn, const std::vector<segment_desc_t> &_segment_desc)
	{
		try {
			traj = new chemfiles::Trajectory(fn);
		}
		catch (...)
		{
			return false;
		}

		try
		{
			ch_frame = traj->read();
		}
		catch(...)
		{
			delete traj;	
			traj = nullptr;
			return false;
		}

		n_atoms = ch_frame.size();
		ch_frame_buffered = true;

//		raw_coords = new rvec[n_atoms];
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
		if (!traj)
			return false;

		try {
			traj->close();
		}
		catch (...)
		{
			return false;
		}

		delete traj;
		traj = nullptr;

		return true;
	}

	bool get_frame(frame_t &frame)
	{
		if (traj == nullptr || traj->done())
			return false;

		frame.clear();
		frame.desc.n_atoms = n_atoms;

		if (!ch_frame_buffered)
		{
			try {
				ch_frame = traj->read();
			}
			catch (...)
			{
				return false;
			}
		}
		else
			ch_frame_buffered = false;

		frame.desc.n_atoms = ch_frame.size();
		frame.desc.step = (int) ch_frame.step();

		if (auto time_opt = ch_frame.get("time"); time_opt)
			frame.desc.time = time_opt->as_double();
		else
			frame.desc.time = 0.0;

		frame.desc.prec = 100;			// for compatibility with XTC precision

		auto matrix = ch_frame.cell().matrix();

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				frame.desc.box[i][j] = matrix[i][j];

		const auto& raw_coords = ch_frame.positions();

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