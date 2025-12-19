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
#include <numeric>
#include <cstring>
#include <iostream>

#include "../common/frame.h"

// *******************************************************************************
class OutputWriter
{
protected:
	int n_atoms{};
	int step{};
	float time{};
	matrix3d_t box{};
	float prec{};
	bool initialized = false;

public:
	OutputWriter()
	{
	}
	virtual ~OutputWriter()
	{
	}

	virtual bool open(const std::string& fn) = 0;
	virtual bool close() = 0;
	virtual bool add_frame(frame_t& frame) = 0;
	virtual bool add_frames(const mdc::query_result& result) = 0;
};

#ifdef USE_XDRFILE
// *******************************************************************************
class CXTCWriter : public OutputWriter
{
	XDRFILE* xtc = nullptr;
	rvec* raw_coords = nullptr;

public:
	CXTCWriter() : OutputWriter()
	{}

	~CXTCWriter()
	{
		if (xtc)
			close();
	}

	virtual bool open(const std::string& fn)
	{
		xtc = xdrfile_open(fn.c_str(), "w");

		if (xtc == nullptr)
			return false;

		initialized = false;

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

	//mkokot
	virtual bool add_frames(const mdc::query_result& result)
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
				if (n_atoms != std::ssize(result.frames[frame_id].coords))
				{
					std::cerr << "when writing to xtc the same number of atoms must be used each time\n";
					exit(1);
				}
			}

			for (int32_t i = 0; i < n_atoms; ++i)
			{
				raw_coords[i][0] = result.frames[frame_id].coords[i].x / 10.0;
				raw_coords[i][1] = result.frames[frame_id].coords[i].y / 10.0;
				raw_coords[i][2] = result.frames[frame_id].coords[i].z / 10.0;
			}

			//all should be the same
			auto& desc = result.frames[frame_id];
			matrix box_float;

			for(int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box_float[i][j] = float(desc.box[i][j]) / 10.0;

			// (desc.prec * 10) must be used since internally we operate in A (not in nm)
			if (write_xtc(xtc, n_atoms, desc.step, desc.time, box_float, raw_coords, desc.prec * 10) != exdrOK)
				return false;
		}
		return true;
	}

	virtual bool add_frame(frame_t& frame)
	{
		if (xtc == nullptr)
			return false;

		if (!initialized)
		{
			n_atoms = 0;
			for (auto& sd : frame.segments)
				n_atoms += (int)sd.coords.size();

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
				raw_coords[raw_i][0] = (float)coords[i].x / frame.desc.prec / 10.0;
				raw_coords[raw_i][1] = (float)coords[i].y / frame.desc.prec / 10.0;
				raw_coords[raw_i][2] = (float)coords[i].z / frame.desc.prec / 10.0;
			}
		}

		matrix box_float;

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				box_float[i][j] = float(frame.desc.box[i][j]) / 10.0;

		// (frame.desc.prec * 10) must be used since internally we operate in A (not in nm)
		if (write_xtc(xtc, n_atoms, frame.desc.step, frame.desc.time, box_float, raw_coords, frame.desc.prec * 10) != exdrOK)
			return false;

		return true;
	}
};

// *******************************************************************************
class CTRRWriter : public OutputWriter
{
	XDRFILE* trr = nullptr;
	rvec* raw_coords = nullptr;

public:
	CTRRWriter() : OutputWriter()
	{
	}

	~CTRRWriter()
	{
		if (trr)
			close();
	}

	virtual bool open(const std::string& fn)
	{
		trr = xdrfile_open(fn.c_str(), "w");

		if (trr == nullptr)
			return false;

		initialized = false;

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

	//mkokot
	virtual bool add_frames(const mdc::query_result& result)
	{
		if (trr == nullptr)
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
				if (n_atoms != std::ssize(result.frames[frame_id].coords))
				{
					std::cerr << "when writing to xtc the same number of atoms must be used each time\n";
					exit(1);
				}
			}

			for (int32_t i = 0; i < n_atoms; ++i)
			{
				raw_coords[i][0] = result.frames[frame_id].coords[i].x / 10.0;
				raw_coords[i][1] = result.frames[frame_id].coords[i].y / 10.0;
				raw_coords[i][2] = result.frames[frame_id].coords[i].z / 10.0;
			}

			//all should be the same
			auto& desc = result.frames[frame_id];
			matrix box_float;

			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					box_float[i][j] = float(desc.box[i][j]) / 10.0;

			if (write_trr(trr, n_atoms, desc.step, desc.time, desc.lambda, box_float, raw_coords, nullptr, nullptr) != exdrOK)
				return false;
		}
		return true;
	}

	virtual bool add_frame(frame_t& frame)
	{
		if (trr == nullptr)
			return false;

		if (!initialized)
		{
			n_atoms = 0;
			for (auto& sd : frame.segments)
				n_atoms += (int)sd.coords.size();

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
				raw_coords[raw_i][0] = (float)coords[i].x / frame.desc.prec / 10.0;
				raw_coords[raw_i][1] = (float)coords[i].y / frame.desc.prec / 10.0;
				raw_coords[raw_i][2] = (float)coords[i].z / frame.desc.prec / 10.0;
			}
		}

		matrix box_float;

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				box_float[i][j] = float(frame.desc.box[i][j]) / 10.0;


		if (write_trr(trr, n_atoms, frame.desc.step, frame.desc.time, frame.desc.lambda, box_float, raw_coords, nullptr, nullptr) != exdrOK)
			return false;

		return true;
	}
};
#endif

#ifdef USE_CHEMFILES
// *******************************************************************************
class CTrajWriter : public OutputWriter
{
	chemfiles::Trajectory* traj = nullptr;
	chemfiles::Frame ch_frame;

public:
	CTrajWriter() : OutputWriter()
	{}

	~CTrajWriter()
	{
		if (traj)
			close();
	}

	virtual bool open(const std::string& fn)
	{
		try {
			traj = new chemfiles::Trajectory(fn, 'w');
		}
		catch (...)
		{
			return false;
		}

		initialized = false;

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

	//mkokot
	virtual bool add_frames(const mdc::query_result& result)
	{
		if (traj == nullptr)
			return false;

		for (uint32_t frame_id = 0 ; frame_id < result.frames.size(); ++frame_id)
		{
			if (!initialized)
			{
				n_atoms = result.frames[frame_id].coords.size();

				ch_frame.resize(n_atoms);

				initialized = true;
			}
			else
			{
				if (n_atoms != std::ssize(result.frames[frame_id].coords))
				{
					std::cerr << "when writing to xtc the same number of atoms must be used each time\n";
					exit(1);
				}
			}

			auto ch_positions = ch_frame.positions();
			for (int32_t i = 0; i < n_atoms; ++i)
			{
				ch_positions[i][0] = result.frames[frame_id].coords[i].x;
				ch_positions[i][1] = result.frames[frame_id].coords[i].y;
				ch_positions[i][2] = result.frames[frame_id].coords[i].z;
			}
			auto& frame = result.frames[frame_id];
			
			ch_frame.set_cell(chemfiles::UnitCell(chemfiles::Matrix3D(
				frame.box[0][0], frame.box[0][1], frame.box[0][2],
				frame.box[1][0], frame.box[1][1], frame.box[1][2],
				frame.box[2][0], frame.box[2][1], frame.box[2][2]
			)));

//			ch_frame.set_index(frame.step);
			ch_frame.set("simulation_step", frame.step);
			ch_frame.set("time", frame.time);
			ch_frame.set("xtc_precision", frame.prec * 10);

			traj->write(ch_frame);
		}

		return true;
	}

	virtual bool add_frame(frame_t& frame)
	{
		if (traj == nullptr)
			return false;

		if (!initialized)
		{
			n_atoms = 0;
			for (auto& sd : frame.segments)
				n_atoms += (int)sd.coords.size();

			frame.desc.n_atoms = n_atoms;

			ch_frame.resize(n_atoms);

			initialized = true;
		}

		uint32_t raw_i = 0;

		for (auto& sd : frame.segments)
		{
			auto& coords = sd.coords;

			auto ch_positions = ch_frame.positions();

			for (uint32_t i = 0; i < sd.coords.size(); ++i, ++raw_i)
			{
				ch_positions[raw_i][0] = (float)coords[i].x / frame.desc.prec;
				ch_positions[raw_i][1] = (float)coords[i].y / frame.desc.prec;
				ch_positions[raw_i][2] = (float)coords[i].z / frame.desc.prec;
			}
		}

		ch_frame.set_cell(chemfiles::UnitCell(chemfiles::Matrix3D(
			frame.desc.box[0][0], frame.desc.box[0][1], frame.desc.box[0][2],
			frame.desc.box[1][0], frame.desc.box[1][1], frame.desc.box[1][2],
			frame.desc.box[2][0], frame.desc.box[2][1], frame.desc.box[2][2]
		)));

//		ch_frame.set_index(frame.desc.step);
		ch_frame.set("simulation_step", frame.desc.step);
		ch_frame.set("time", frame.desc.time);
		ch_frame.set("xtc_precision", frame.desc.prec * 10);

		traj->write(ch_frame);

		return true;
	}
};
#endif