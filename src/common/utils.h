#pragma once

#include <string>
#include <algorithm>
#include <cmath>

template<typename D, typename T>
D iround(T x)
{
	if (x < 0)
		return D(x - 0.5);
	else
		return D(x + 0.5);
}

template<typename D, typename T>
D uround(T x)
{
//	return (D)(x + 0.5);
	return static_cast<D>(round(x));
}

template<typename D, typename T>
D uround_div(T x, T y)
{
	return static_cast<D>((x + y / 2) / y);
}

template<typename T>
T iabs(T x)
{
	if (x < 0)
		return -x;
	else
		return x;
}

template<typename T>
T idiv(T x, T y)
{
	bool negative = false;

	if (y < 0)
	{
		y = -y;
		negative = true;
	}

	if (x < 0)
	{
		x = -x;
		negative = !negative;
	}

	T r = x / y;
	
	if (x - r * y >= y / 2)
		r++;

	return negative ? -r : r;
}

enum class traj_file_format_t {unknown, trr, xtc, tng, mdc, nc, gro, dcd};

constexpr std::vector<std::string> list_of_supported_input_formats_for_compression()
{
	return { "trr", "xtc", "tng","nc", "gro", "dcd" };
}

constexpr std::vector<std::string> list_of_supported_output_formats_for_decompression()
{
	return { "trr", "xtc", "nc", "gro", "dcd" };
}

constexpr traj_file_format_t get_traj_format(const std::string& fn, bool for_output)
{
	if (fn.empty())
		return traj_file_format_t::unknown;

	int i = int(fn.size()) - 1;
	for (; i >= 0; --i)
		if (fn[i] == '.')
			break;

	if (i < 0)
		return traj_file_format_t::unknown;

	auto ext = fn.substr(i + 1);

	if (ext == "trr" || ext == "TRR")
		return traj_file_format_t::trr;
	if (ext == "xtc" || ext == "XTC")
		return traj_file_format_t::xtc;
	if (ext == "nc" || ext == "NC")
		return traj_file_format_t::nc;
	if (ext == "mdc" || ext == "MDC")
		return traj_file_format_t::mdc;
	if (ext == "gro" || ext == "GRO")
		return traj_file_format_t::gro;
	if (ext == "dcd" || ext == "DCD")
		return traj_file_format_t::dcd;

	if(for_output)
		return traj_file_format_t::unknown;

	if (ext == "tng" || ext == "TNG")
		return traj_file_format_t::tng;

	return traj_file_format_t::unknown;
}