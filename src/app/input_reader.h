#pragma once

#include "defs.h"

#ifdef USE_XDRFILE
#include "../../libs/libxdrfile/xdrfile_xtc.h"
#include "../../libs/libxdrfile/xdrfile_trr.h"
#endif

#if defined(ONLY_XTC_READER) && defined(USE_CHEMFILES)
#undef USE_CHEMFILES
#endif

#ifdef USE_CHEMFILES
#include <chemfiles.hpp>
#endif

#include <string>
#include <span>

#include "../common/utils.h"
#include "../common/frame.h"

// *******************************************************************************
class InputReader
{
protected:
	int n_atoms{};
	int step{};
	float time{};
	matrix3d_t box{};
	float prec{};
	bool initialized = false;
	uint32_t resolution = 1000; // in fm

	std::vector<segment_desc_t> segment_desc;
	std::vector<icoord_t> coords;

public:
	InputReader(uint32_t resolution) :
		resolution(resolution)
	{}
	virtual ~InputReader()
	{}

	virtual bool open(const std::string& fn, const std::vector<segment_desc_t>& _segment_desc) = 0;
	virtual bool close() = 0;
	virtual bool get_frame(frame_t& frame) = 0;

	void adjust_precision()
	{
		if(resolution != 0)
			prec = (float)100000.0 / resolution;
	}

	uint32_t get_resolution() const
	{
		return resolution;
	}

	template<typename T>
	void round_coords(const std::span<T> raw_coords, float mult = 1.0)
	{
		coords.resize(raw_coords.size());

//		auto loc_prec = prec * mult;

		if (mult != 1.0)
			for (auto i = 0; i < std::ssize(raw_coords); ++i)
			{
				//			coords[i].x = iround<int>(raw_coords[i][0] * loc_prec);
				//			coords[i].y = iround<int>(raw_coords[i][1] * loc_prec);
				//			coords[i].z = iround<int>(raw_coords[i][2] * loc_prec);
				coords[i].x = iround<int>((raw_coords[i][0] * mult) * prec);
				coords[i].y = iround<int>((raw_coords[i][1] * mult) * prec);
				coords[i].z = iround<int>((raw_coords[i][2] * mult) * prec);
			}
		else
			for (auto i = 0; i < std::ssize(raw_coords); ++i)
			{
				coords[i].x = iround<int>(raw_coords[i][0] * prec);
				coords[i].y = iround<int>(raw_coords[i][1] * prec);
				coords[i].z = iround<int>(raw_coords[i][2] * prec);
			}
	}
};

#ifdef USE_XDRFILE
// *******************************************************************************
class CXTCReader : public InputReader
{
	XDRFILE* xtc = nullptr;
	rvec* raw_coords = nullptr;

public:
	CXTCReader(uint32_t resolution) : InputReader(resolution)
	{
	}

	virtual ~CXTCReader()
	{
		if (xtc)
			close();
	}

	virtual bool open(const std::string& fn, const std::vector<segment_desc_t>& _segment_desc)
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

		if (n_atoms_in_desc != (uint32_t)n_atoms)
			return false;

		return true;
	}

	virtual bool close()
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

	virtual bool get_frame(frame_t& frame)
	{
		if (xtc == nullptr)
			return false;

		frame.clear();
		frame.desc.n_atoms = n_atoms;

		matrix box_float;

		if (read_xtc(xtc, frame.desc.n_atoms, &frame.desc.step, &frame.desc.time, box_float, raw_coords, &frame.desc.prec) != exdrOK)
			return false;

		for(int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				frame.desc.box[i][j] = box_float[i][j] * 10.0;

		if (!initialized)
		{
			step = frame.desc.step;
			time = frame.desc.time;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box[i][j] = frame.desc.box[i][j];
			prec = frame.desc.prec / 10.0;					// all coordinates are processed in A (for compatibility with chemfiles)

			initialized = true;

			adjust_precision();
		}

		frame.desc.prec = prec;								// adjusted precision will be used consistently
		
		round_coords(std::span(raw_coords, n_atoms), 10);

		int raw_i = 0;

		for (uint32_t i = 0; i < (uint32_t)segment_desc.size(); ++i)
		{
			uint32_t n_seg_atoms = segment_desc[i].size;

			// Do not load segments marked as NONE
			if (segment_desc[i].type == segment_type_t::none)
			{
				raw_i += n_seg_atoms;
				continue;
			}

			segment_t segment(segment_desc[i].type, segment_desc[i].name, n_seg_atoms, std::vector<icoord_t>());

			auto& s_coords = segment.coords;

			s_coords.assign(coords.begin() + raw_i, coords.begin() + raw_i + n_seg_atoms);
			raw_i += n_seg_atoms;
			/*			s_coords.clear();
						s_coords.reserve(n_seg_atoms);

						for (uint32_t j = 0; j < n_seg_atoms; ++j, ++raw_i)
							s_coords.emplace_back(coords[raw_i]);
							*/

			frame.emplace(std::move(segment));
		}

		return true;
	}
};

// *******************************************************************************
class CTRRReader : public InputReader
{
	XDRFILE* trr = nullptr;
	rvec* raw_coords = nullptr;

public:
	CTRRReader(uint32_t resolution) : InputReader(resolution)
	{
	}

	virtual ~CTRRReader()
	{
		if (trr)
			close();
	}

	virtual bool open(const std::string& fn, const std::vector<segment_desc_t>& _segment_desc)
	{
		if (read_trr_natoms((char*)fn.c_str(), &n_atoms) != exdrOK)
			return false;

		trr = xdrfile_open(fn.c_str(), "r");

		if (trr == nullptr)
			return false;

		raw_coords = new rvec[n_atoms];
		initialized = false;

		segment_desc = _segment_desc;

		uint32_t n_atoms_in_desc = 0;
		for (auto& sd : segment_desc)
			n_atoms_in_desc += sd.size;

		if (n_atoms_in_desc != (uint32_t)n_atoms)
			return false;

		return true;
	}

	virtual bool close()
	{
		if (!trr)
			return false;

		if (xdrfile_close(trr) != exdrOK)
			return false;

		trr = nullptr;

		delete[] raw_coords;
		raw_coords = nullptr;

		return true;
	}

	virtual bool get_frame(frame_t& frame)
	{
		if (trr == nullptr)
			return false;

		frame.clear();
		frame.desc.n_atoms = n_atoms;

		matrix box_float;

		if (read_trr(trr, frame.desc.n_atoms, &frame.desc.step, &frame.desc.time, &frame.desc.lambda, box_float, 
			raw_coords, nullptr, nullptr, &frame.desc.has_prop) != exdrOK)
			return false;

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				frame.desc.box[i][j] = box_float[i][j] * 10.0;

		if (!initialized)
		{
			step = frame.desc.step;
			time = frame.desc.time;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box[i][j] = frame.desc.box[i][j];
			prec = 1e5;

			adjust_precision();

			initialized = true;
		}

		frame.desc.prec = prec;								// adjusted precision will be used consistently

		round_coords(std::span(raw_coords, n_atoms), 10);

		int raw_i = 0;

		for (uint32_t i = 0; i < (uint32_t)segment_desc.size(); ++i)
		{
			uint32_t n_seg_atoms = segment_desc[i].size;

			// Do not load segments marked as NONE
			if (segment_desc[i].type == segment_type_t::none)
			{
				raw_i += n_seg_atoms;
				continue;
			}

			segment_t segment(segment_desc[i].type, segment_desc[i].name, n_seg_atoms, std::vector<icoord_t>());

			auto& s_coords = segment.coords;

			s_coords.assign(coords.begin() + raw_i, coords.begin() + raw_i + n_seg_atoms);
			raw_i += n_seg_atoms;
			/*			s_coords.clear();
						s_coords.reserve(n_seg_atoms);

						for (uint32_t j = 0; j < n_seg_atoms; ++j, ++raw_i)
							s_coords.emplace_back(coords[raw_i]);
							*/

			frame.emplace(std::move(segment));
		}

		return true;
	}
};
#endif

#ifdef USE_CHEMFILES
// *******************************************************************************
class CTrajReader : public InputReader
{
	chemfiles::Trajectory* traj = nullptr;
	chemfiles::Frame ch_frame;

	bool ch_frame_buffered;

public:
	CTrajReader(uint32_t no_fpd_A) : InputReader(no_fpd_A)
	{
	}

	virtual ~CTrajReader()
	{
		if (traj)
			close();
	}

	virtual bool open(const std::string& fn, const std::vector<segment_desc_t>& _segment_desc)
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
		catch (...)
		{
			delete traj;
			traj = nullptr;
			return false;
		}

		n_atoms = ch_frame.size();
		ch_frame_buffered = true;

		initialized = false;

		segment_desc = _segment_desc;

		uint32_t n_atoms_in_desc = 0;
		for (auto& sd : segment_desc)
			n_atoms_in_desc += sd.size;

		if (n_atoms_in_desc != (uint32_t)n_atoms)
			return false;

		return true;
	}

	virtual bool close()
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

	virtual bool get_frame(frame_t& frame)
	{
		if (traj == nullptr || traj->done())
			return false;

		frame.clear();
		frame.desc.n_atoms = n_atoms;

		if (!ch_frame_buffered)
		{
			try {
				ch_frame = std::move(traj->read());
			}
			catch (...)
			{
				return false;
			}
		}
		else
			ch_frame_buffered = false;

		frame.desc.n_atoms = ch_frame.size();

		if (auto time_opt = ch_frame.get("time"); time_opt)
			frame.desc.time = time_opt->as_double();
		else
			frame.desc.time = 0.0;

//		frame.desc.prec = 100;			// for compatibility with XTC precision
		if (auto xtc_precision = ch_frame.get("xtc_precision"); xtc_precision)
			frame.desc.prec = xtc_precision->as_double() / 10.0;
		else
			frame.desc.prec = 100;

		if (auto simulation_step = ch_frame.get("simulation_step"); simulation_step)
			frame.desc.step = (int) simulation_step->as_double();
		else
			frame.desc.step = (int) ch_frame.index();

//		frame.desc.step = (int)ch_frame.index();

		auto matrix = ch_frame.cell().matrix();

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				frame.desc.box[i][j] = matrix[i][j];

		if (!initialized)
		{
			step = frame.desc.step;
			time = frame.desc.time;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box[i][j] = frame.desc.box[i][j];
			prec = frame.desc.prec;

			adjust_precision();

			initialized = true;
		}

		frame.desc.prec = prec;								// adjusted precision will be used consistently

		const auto& ch_coords = ch_frame.positions();

		round_coords(std::span(ch_coords.begin(), ch_coords.end()));

		if (!initialized)
		{
			step = frame.desc.step;
			time = frame.desc.time;
			box = frame.desc.box;
			prec = frame.desc.prec;

			initialized = true;
		}

		int raw_i = 0;

		for (uint32_t i = 0; i < (uint32_t)segment_desc.size(); ++i)
		{
			uint32_t n_seg_atoms = segment_desc[i].size;

			// Do not load segments marked as NONE
			if (segment_desc[i].type == segment_type_t::none)
			{
				raw_i += n_seg_atoms;
				continue;
			}

			segment_t segment(segment_desc[i].type, segment_desc[i].name, n_seg_atoms, std::vector<icoord_t>());

			auto& s_coords = segment.coords;

			s_coords.assign(coords.begin() + raw_i, coords.begin() + raw_i + n_seg_atoms);
			raw_i += n_seg_atoms;
/*			s_coords.clear();
			s_coords.reserve(n_seg_atoms);

			for (uint32_t j = 0; j < n_seg_atoms; ++j, ++raw_i)
				s_coords.emplace_back(coords[raw_i]);
				*/
			frame.emplace(std::move(segment));
		}

		return true;
	}
};
#endif
