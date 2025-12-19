#pragma once

#include <vector>
#include <array>
#include <memory>

#include "frame.h"
#include "serializer.h"
#include "engine.h"
#include "../app/params.h"

// *******************************************************************************
class CSegmentModel
{
public:
	CSegmentModel() {}
	virtual ~CSegmentModel() {}
};

// *******************************************************************************
class CSegmentModelOther : public CSegmentModel
{
public:
	CSegmentModelOther() : CSegmentModel() {}
	virtual ~CSegmentModelOther() {}
};

// *******************************************************************************
class CSegmentModelWater : public CSegmentModel
{
public:
	uint64_t dist2_OH{};
	uint64_t dist2_HH{};
	uint8_t atom_order{};	// 0 - OHH, 1 - HOH, 2 - HHO

public:
	CSegmentModelWater() : CSegmentModel() {}
	virtual ~CSegmentModelWater() {}
};

// *******************************************************************************
class CSegmentModelMolecule : public CSegmentModel
{
public:
	struct ref_atoms_t
	{
		std::array<uint32_t, 3> atom_ids{};
		std::array<uint32_t, 3> distances{};
		bool is_left_handed = false;

		ref_atoms_t() = default;

		ref_atoms_t(uint32_t id0, uint32_t id1, uint32_t id2, uint32_t d0, uint32_t d1, uint32_t d2, bool _is_left_handed = false)
		{
			atom_ids[0] = id0;
			atom_ids[1] = id1;
			atom_ids[2] = id2;
			distances[0] = d0;
			distances[1] = d1;
			distances[2] = d2;
			is_left_handed = _is_left_handed;
		}
	};

	std::vector<ref_atoms_t> ref_atoms;

public:
	CSegmentModelMolecule() : CSegmentModel() {}
	virtual ~CSegmentModelMolecule() {}
};

// *******************************************************************************
class CSegmentHistory
{
public:
	uint32_t max_history_size{};
	std::vector<segment_t> history;

	CSegmentHistory(const CSegmentHistory&) = default;
	CSegmentHistory(CSegmentHistory&&) noexcept = default;

	CSegmentHistory(uint32_t max_history_size = 1) :
		max_history_size(max_history_size)
	{}

	CSegmentHistory& operator=(const CSegmentHistory&) = default;
	CSegmentHistory& operator=(CSegmentHistory&&) noexcept = default;

	void add(const segment_t& segment)
	{
		if (history.size() >= max_history_size)
			history.pop_back();
		history.emplace(history.begin(), segment);
	}

	void add(segment_t&& segment)
	{
		if (history.size() >= max_history_size)
			history.pop_back();
		history.emplace(history.begin(), std::move(segment));
	}

	void clear()
	{
		history.clear();
	}

	bool empty() const
	{
		return history.empty();
	}

	size_t size() const
	{
		return history.size();
	}

	size_t size(uint32_t i) const
	{
		return history[i].coords.size();
	}	

	void get(uint32_t i, icoord_t &a) const
	{
		a = history[0].coords[i];
	}

	void get(uint32_t i, icoord_t& a, icoord_t &b) const
	{
		a = history[0].coords[i];
		b = history[1].coords[i];
	}

	void get(uint32_t i, icoord_t& a, icoord_t &b, icoord_t &c) const
	{
		a = history[0].coords[i];
		b = history[1].coords[i];
		c = history[2].coords[i];
	}

	double dist(uint32_t i, uint32_t j) const
	{
		auto &a = history[0].coords[i];
		auto &b = history[0].coords[j];
		auto dx = a.x - b.x;
		auto dy = a.y - b.y;
		auto dz = a.z - b.z;

		return sqrt(dx * dx + dy * dy + dz * dz);
	}

	double dist2(uint32_t i, uint32_t j) const
	{
		auto &a = history[0].coords[i];
		auto &b = history[0].coords[j];
		auto dx = a.x - b.x;
		auto dy = a.y - b.y;
		auto dz = a.z - b.z;

		return dx * dx + dy * dy + dz * dz;
	}
};

// *******************************************************************************
class CMDCAbstract
{
public:
	struct packed_frames_t
	{
		packed_t desc;
		std::vector<packed_t> packed;

		packed_frames_t() = default;
		packed_frames_t(const packed_frames_t&) = default;
		packed_frames_t(packed_frames_t&&) = default;
		packed_frames_t& operator=(const packed_frames_t&) = default;
		packed_frames_t& operator=(packed_frames_t&&) = default;

		bool empty() const
		{
			return packed.empty();
		}
	};

protected:
	enum class frame_type_t {anchor = 0, delta = 1};

	uint32_t no_delta_frames{0};
	uint32_t max_history_size{0};

	compression_features_t compression_features;

	Serializer serializer;

	std::vector<uint32_t> anchor_ids;

	std::vector<segment_desc_t> subsegment_desc;
	
	std::vector<std::shared_ptr<CSegmentModel>> segment_models;
	std::vector<CSegmentHistory> segment_histories;

	std::vector<frame_type_t> frame_types;

public:
	CMDCAbstract(uint32_t no_delta_frames = 0, uint32_t max_history_size = 0, compression_features_t compression_features = compression_features_t()) :
		no_delta_frames(no_delta_frames),
		max_history_size(max_history_size),
		compression_features(compression_features)
	{}

	virtual ~CMDCAbstract()
	{}
};
