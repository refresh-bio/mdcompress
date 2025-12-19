#pragma once

#include <string>
#include <cinttypes>
#include <thread>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>
#include <cassert>

enum class preset { default_, archive, trajectory, frames};

inline std::optional<preset> preset_from_string(const std::string& s)
{
	if (s == "default")
		return preset::default_;
	if (s == "archive")
		return preset::archive;
	if (s == "trajectory")
		return preset::trajectory;
	if (s == "frames")
		return preset::frames;

	return std::nullopt;
};

inline std::string to_string(preset p)
{
	switch (p)
	{
		case preset::default_: return "default";
		case preset::archive: return "archive";
		case preset::trajectory: return "trajectory";
		case preset::frames: return "frames";
	}
	return "unknown";
}


// *******************************************************************************************
template<typename T> class b_value
{
	T value;
	T min_value;
	T max_value;

	std::string trim_trailing_zeros(const std::string& str) const
	{
		std::string result = str;
		if (result.find('.') != std::string::npos)
		{
			result.erase(result.find_last_not_of('0') + 1, std::string::npos);
			if (result.back() == '.')
				result.pop_back();
		}
		return result;
	}

public:
	b_value(const T _value, const T _min_value, const T _max_value) : value(_value), min_value(_min_value), max_value(_max_value) {};

	bool assign(const T _value)
	{
		value = std::clamp(_value, min_value, max_value);

		return value == _value;
	}

	std::string info() const
	{
		std::string valueStr = to_string(value);
		std::string minStr = to_string(min_value);
		std::string maxStr = to_string(max_value);

		if (std::is_floating_point<T>::value)
		{
			valueStr = trim_trailing_zeros(valueStr);
			minStr = trim_trailing_zeros(minStr);
			maxStr = trim_trailing_zeros(maxStr);
		}

		/*		return "(default: " + valueStr +
					"; min: " + minStr +
					"; max: " + maxStr + ")";*/

		return "(default: " + valueStr +
			"; range " + minStr +
			"-" + maxStr + ")";
	}

	std::string range(int offset = 0) const
	{
		return "[" + std::to_string(min_value + offset) + "-" + std::to_string(max_value + offset) + "]";
	}

	std::string info_no_range() const
	{
		std::string valueStr = to_string(value);
		std::string minStr = to_string(min_value);
		std::string maxStr = to_string(max_value);

		if (std::is_floating_point<T>::value)
		{
			valueStr = trim_trailing_zeros(valueStr);
			minStr = trim_trailing_zeros(minStr);
			maxStr = trim_trailing_zeros(maxStr);
		}

		/*		return "(default: " + valueStr +
					"; min: " + minStr +
					"; max: " + maxStr + ")";*/

		return "(default: " + valueStr + ")";
	}

	T operator()() const
	{
		return value;
	}

	T min_val() const
	{
		return min_value;
	}

	T max_val() const
	{
		return max_value;
	}
};

// *******************************************************************************************
//value that will be always divisible by DIVISOR
//this class will auto adjust it
template<typename T, unsigned DIVISOR>
class b_value_divisible
{
	static_assert(DIVISOR > 0);
	static_assert(std::is_integral_v<T>);
	b_value<T> value;
	static T make_divisible(T v)
	{
		return (v + DIVISOR - 1) / DIVISOR * DIVISOR;
	}
public:
	b_value_divisible(const T _value, const T _min_value, const T _max_value) : value(make_divisible(_value), make_divisible(_min_value), make_divisible(_max_value))
	{

	}

	bool assign(const T _value)
	{
		return value.assign(make_divisible(_value));
	}

	T operator()() const
	{
		return value();
	}

	T min_val() const
	{
		return value.min_val();
	}

	T max_val() const
	{
		return max_val();
	}
	std::string range(int offset = 0) const
	{
		return value.range(offset);
	}
	static consteval unsigned divisor()
	{
		return DIVISOR;
	}
};

struct preset_values_t
{
	//only get via static methods
private:
	preset_values_t(uint32_t anchor_separation,
		uint32_t subsegment_size) :
	anchor_separation(anchor_separation, 0, 999'999),
	subsegment_size(subsegment_size, 0, 1'000'000)  // 0 means do not use subsegments
	{}
	static preset_values_t get_default()
	{
		return { 19, 100 };
	}
	static preset_values_t get_archive()
	{
		return { 99, 1000 };
	}
	static preset_values_t get_trajectory()
	{
		return { 99, 100 };
	}
	static preset_values_t get_frames()
	{
		return { 9, 1000 };
	}
public:
	b_value<uint32_t> anchor_separation;

	//get closest >= value divisible by 3, required for water, does not matter for others, so for simplicity we do it for all
	b_value_divisible<uint32_t, 3> subsegment_size;

	static preset_values_t get(preset p)
	{
		switch (p)
		{
			case preset::default_: return get_default();
			case preset::archive: return get_archive();
			case preset::trajectory: return get_trajectory();
			case preset::frames: return get_frames();
		}
		assert(false);
	}
};


// *******************************************************************************************
struct compression_features_t
{
	enum class water_z1_t { none = 0, predicted = 1 };
	enum class water_x2_t { none = 0, predicted = 1 };
	enum class water_y2_t { none = 0, predicted = 1 };
	enum class water_z2_t { none = 0, predicted = 1 };
	enum class molecule_t { moving_1vector = 0, moving_3vectors = 1, tetrahedron = 2, translate_and_rotate = 3};
	enum class molecule_orientation_t { none = 0, predicted = 1};
	enum class rc_uint32_t { prefix_fast = 0, prefix_exact = 1, suffix_exact = 2};


//	enum class 

	water_z1_t water_z1{ water_z1_t::none };
	water_x2_t water_x2{ water_x2_t::none };
	water_y2_t water_y2{ water_y2_t::none };
	water_z2_t water_z2{ water_z2_t::none };
	molecule_t molecule{ molecule_t::moving_1vector };
	molecule_orientation_t molecule_orientation{ molecule_orientation_t::none };
	rc_uint32_t rc_uint32{ rc_uint32_t::prefix_fast };

	compression_features_t() = default;
	compression_features_t(water_z1_t z1, water_x2_t x2, water_y2_t y2, water_z2_t z2, molecule_t m, molecule_orientation_t mo, rc_uint32_t r) :
		water_z1(z1), water_x2(x2), water_y2(y2), water_z2(z2), molecule(m), molecule_orientation(mo), rc_uint32(r)
	{}
};

// *******************************************************************************************
struct CParamsMakeDesc
{
	std::string input_path;
	std::string output_path;
};

// *******************************************************************************************
struct CParams
{
private:
	static uint32_t def_no_threads()
	{
		const auto hc = std::thread::hardware_concurrency();
		const auto res = std::min(8u, hc);
		return std::max(res, 1u);
	}
public:
	enum class mode_t { unknown, compress, decompress, select, make_desc, convert, info, to_fp32 };

	CParamsMakeDesc params_make_desc;

	using cf = compression_features_t;

	std::vector<compression_features_t> compression_features = { 
		compression_features_t(),				// 0
		compression_features_t{cf::water_z1_t::predicted, cf::water_x2_t::none, cf::water_y2_t::none, cf::water_z2_t::predicted, cf::molecule_t::moving_1vector, cf::molecule_orientation_t::none, cf::rc_uint32_t::prefix_fast },					// 1
		compression_features_t{cf::water_z1_t::predicted, cf::water_x2_t::none, cf::water_y2_t::none, cf::water_z2_t::predicted, cf::molecule_t::moving_1vector, cf::molecule_orientation_t::none, cf::rc_uint32_t::prefix_exact },					// 2
		compression_features_t{cf::water_z1_t::predicted, cf::water_x2_t::none, cf::water_y2_t::predicted, cf::water_z2_t::predicted, cf::molecule_t::moving_1vector, cf::molecule_orientation_t::none, cf::rc_uint32_t::prefix_exact },				// 3
		compression_features_t{cf::water_z1_t::predicted, cf::water_x2_t::none, cf::water_y2_t::predicted, cf::water_z2_t::predicted, cf::molecule_t::translate_and_rotate, cf::molecule_orientation_t::none, cf::rc_uint32_t::prefix_exact },	// 4
		compression_features_t{cf::water_z1_t::predicted, cf::water_x2_t::predicted, cf::water_y2_t::predicted, cf::water_z2_t::predicted, cf::molecule_t::translate_and_rotate, cf::molecule_orientation_t::predicted, cf::rc_uint32_t::suffix_exact }	// 5
	};

	mode_t mode{ mode_t::unknown };
	std::string input_fn{};
	std::string output_fn{};
	std::string description_fn{};
	std::string topology_fn{}; //optional topology file, if present will be stored in the archive, if describtion not provided describtion will be inffered from topology (when whe know how)
	bool only_mol = false;
	b_value<uint32_t> no_threads{ def_no_threads(), 1, 16 };
	int32_t frame_id{ -1 };
	int32_t range_id_from{ -1 };
	int32_t range_id_to{ -1 };
	bool native_lib_mode{ false };
	b_value<uint32_t> resolution{ 1'000, 100, 1'000'000 };
	bool resolution_set_by_user{ false };

	b_value<uint32_t> compression_level{ 2, 1, 5 };
	//we have a couple of preset values, but each may be overriden
	preset_values_t preset_values{ preset_values_t::get(preset::default_) };

	b_value<uint32_t> stride{ 1, 1, 1'000'000'000 };
	b_value<uint32_t> n_frames_for_model{ 50, 1, 1'000 };
	std::vector<std::string> subsegment_ids;
	b_value<uint32_t> max_dist_in_model{ 100, 10, 250 };		// stored at 1 byte, so cannot be larger than 255
	b_value<uint32_t> max_history_size{ 1, 1, 3 };			// give user option to change
	std::vector<std::string> segment_ids;
	std::vector<uint32_t> atom_ids;

	bool full_info = false; //for info mode, prinst also subsegments

	// *******************************************************************************
	std::vector<std::string> split_string(const std::string& s, char sep)
	{
		auto p = s.begin();
		std::vector<std::string> components;

		while (true)
		{
			auto q = find(p, s.end(), sep);
			components.push_back(std::string(p, q));
			if (q == s.end())
				break;

			p = q + 1;
		}

		return components;
	}
};
