#pragma once

#include <cinttypes>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <bit>

#include "../libs/libxdrfile/xdrfile_xtc.h"
#include "../mdc_lib/mdc_reader.h"
#include "utils.h"

//using namespace std;

// *******************************************************************************
template<typename T>
struct coord_t
{
	T x, y, z;

	coord_t() :
		x(0), y(0), z(0)
	{}

	coord_t(T x, T y, T z) :
		x(x), y(y), z(z)
	{}

	coord_t(const coord_t&) = default;
	coord_t(coord_t&&) = default;

	coord_t& operator=(const coord_t&) = default;
	coord_t& operator=(coord_t&) = default;

	bool operator<(const coord_t &rhs)
	{
		if (x != rhs.x)
			return x < rhs.x;
		if (y != rhs.y)
			return y < rhs.y;
		return z < rhs.z;
	}
};

using icoord_t = coord_t<int32_t>;
using fcoord_t = coord_t<float>;
using dcoord_t = coord_t<double>;

// *******************************************************************************
using matrix3d_t = std::array<std::array<double, 3>, 3>;

// *******************************************************************************
template<typename T>
inline bool operator<(const coord_t<T>& a, const coord_t<T>& b)
{
	if (a.x != b.x)
		return a.x < b.x;
	if (a.y != b.y)
		return a.y < b.y;
	return a.z < b.z;
}

// *******************************************************************************
inline uint64_t dist2_3d(const icoord_t& a, const icoord_t& b)
{
	return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
}

// *******************************************************************************
// Function sqrt()
template <typename T>
T isqrt(T x)
{
	if (x == 0)
		return 0;

	int bw = std::bit_width(x);

	T s = T(1) << ((bw + 1) / 2);

	// A few Newton iterations
	// s_{k+1} = floor((s_k + X/s_k)/2)
	for (int i = 0; i < 4; i++) {
		T d = x / s;
		T s_next = (s + d) >> 1;
		s = s_next;
	}

	// For sure round up or down
	while (s * s > x)
		s--;
	while ((s + 1) * (s + 1) <= x)
		s++;

	T s_sq = s * s;

	s += (x - s_sq) > (s_sq + 2 * s + 1 - x);

	return s;
}

// *******************************************************************************
inline uint64_t dist_3d(const icoord_t& a, const icoord_t& b)
{
	return isqrt(dist2_3d(a, b));
}

// *******************************************************************************
using segment_type_t = mdc::segment_type;

using packed_t = std::vector<uint8_t>;

// *******************************************************************************
using segment_desc_t = mdc::segment_desc_t;

// *******************************************************************************
struct segment_t
{
	segment_type_t type;
	std::string name;
	uint32_t n_atoms;

	packed_t packed;
	std::vector<icoord_t> coords;

	segment_t() :
		type(segment_type_t::unknown),
		name(""),
		n_atoms(0)
	{}

	segment_t(const segment_type_t type, const std::string &name, const uint32_t n_atoms, const packed_t &packed) :
		type(type),
		name(name),
		n_atoms(n_atoms),
		packed(packed)
	{}

	segment_t(const segment_type_t type, const std::string& name, const uint32_t n_atoms, const std::vector<icoord_t> &coords) :
		type(type),
		name(name),
		n_atoms(n_atoms),
		coords(coords)
	{}

	segment_t(const segment_t&) = default;
	segment_t(segment_t&&) noexcept  = default;

	segment_t& operator=(const segment_t&) = default;
	segment_t& operator=(segment_t&&) noexcept = default;
};

// *******************************************************************************
struct frame_desc_t
{
	int n_atoms{};
	int step{};
	float time{};
	float lambda{};
	matrix3d_t box{};
	float prec{};
	int has_prop{};

	void clear()
	{
		n_atoms = 0;
		step = 0;	
		time = 0.0f;
		for(int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				box[i][j] = 0.0f;
		prec = 0;

		lambda = 0;
		has_prop = 0;
	}
}; 

struct frame_t
{
	bool compressed;

	std::vector<segment_t> segments; //mkokot_TODO: should rename to subsegment? probably yes because this type is not exposed to the user so we can name it in internal/non-user nomenclature (user should not be notified of subsegments)
	frame_desc_t desc;

	frame_t() :
		compressed(false)
	{}

	frame_t(const frame_t&) = default;
	frame_t(frame_t&&) = default;

	frame_t& operator=(const frame_t&) = default;
	frame_t& operator=(frame_t&&) = default;

	void append(const segment_t& segment)
	{
		segments.emplace_back(segment);
	}

	void emplace(segment_t&& segment)
	{
		segments.emplace_back(std::move(segment));
	}

	void add(const uint32_t segment_id, const segment_t& segment)
	{
		if (segment_id >= segments.size())
			segments.resize(segment_id+1);

		segments[segment_id] = segment;
	}

	void clear()
	{
		compressed = false;
		segments.clear();
	}

	void shrink_to_fit()
	{
		segments.shrink_to_fit();
	}
};

