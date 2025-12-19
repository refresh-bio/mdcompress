#pragma once

#include <cmath>
#include "../libs/refresh/fp_float/lib/fp_float.h"

// *****************************************************************
struct Vec3Dfp
{
	using my_float = refresh::fp_float<int64_t, 16>;

	my_float x, y, z, t;

	Vec3Dfp() : x(0), y(0), z(0)
	{}

	Vec3Dfp(my_float x, my_float y, my_float z) : x(x), y(y), z(z)
	{}

	Vec3Dfp operator-() const
	{
		return Vec3Dfp(-x, -y, -z);
	}

	Vec3Dfp& operator+=(const Vec3Dfp& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;

		return *this;
	}

	bool operator<(const Vec3Dfp& rhs) const
	{
		if (x != rhs.x)
			return x < rhs.x;
		if (y != rhs.y)
			return y < rhs.y;
		return z < rhs.z;
	}

	bool is_wrong() const
	{
		return x.is_wrong() || y.is_wrong() || z.is_wrong();
	}

	friend my_float vector_length(const Vec3Dfp& a);
	friend Vec3Dfp operator-(const Vec3Dfp& a, const Vec3Dfp& b);
	friend Vec3Dfp sub(const Vec3Dfp& a, const Vec3Dfp& b);
	friend Vec3Dfp operator+(const Vec3Dfp& a, const Vec3Dfp& b);
	friend Vec3Dfp operator*(my_float d, const Vec3Dfp& a);
	friend Vec3Dfp operator/(const Vec3Dfp& a, my_float d);
	friend Vec3Dfp normalize_vector(const Vec3Dfp& a);
	friend my_float dot_product(const Vec3Dfp& a, const Vec3Dfp& b);
	friend Vec3Dfp perp_vector(const Vec3Dfp& a, const Vec3Dfp& b);
	friend Vec3Dfp cross_product(const Vec3Dfp& a, const Vec3Dfp& b);
};

inline Vec3Dfp::my_float vector_length(const Vec3Dfp & a)
{
	return refresh::sqrt(pow2(a.x) + pow2(a.y) + pow2(a.z));
//	return refresh::approx_sqrt(pow2(a.x) + pow2(a.y) + pow2(a.z), 10);
}

inline Vec3Dfp::my_float vector_length2(const Vec3Dfp& a)
{
	return pow2(a.x) + pow2(a.y) + pow2(a.z);
}

inline Vec3Dfp operator-(const Vec3Dfp& a, const Vec3Dfp& b)
{
	return Vec3Dfp(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vec3Dfp sub(const Vec3Dfp& a, const Vec3Dfp& b)
{
	return Vec3Dfp(sub(a.x, b.x), sub(a.y, b.y), sub(a.z, b.z));
}

inline Vec3Dfp operator+(const Vec3Dfp& a, const Vec3Dfp& b)
{
	return Vec3Dfp(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3Dfp operator*(Vec3Dfp::my_float d, const Vec3Dfp& a)
{
	return Vec3Dfp(d * a.x, d * a.y, d * a.z);
}

inline Vec3Dfp operator/(const Vec3Dfp& a, Vec3Dfp::my_float d)
{
	return Vec3Dfp(a.x / d, a.y / d, a.z / d);
}

inline Vec3Dfp normalize_vector(const Vec3Dfp& a)
{
	Vec3Dfp::my_float len = vector_length(a);

	return Vec3Dfp(a.x / len, a.y / len, a.z / len);
}

inline Vec3Dfp::my_float dot_product(const Vec3Dfp& a, const Vec3Dfp& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3Dfp perp_vector(const Vec3Dfp& a, const Vec3Dfp& b)
{
	Vec3Dfp::my_float dp = dot_product(a, b);

	return Vec3Dfp(a.x - b.x * dp, a.y - b.y * dp, a.z - b.z * dp);
}

inline Vec3Dfp cross_product(const Vec3Dfp& a, const Vec3Dfp& b)
{
	return Vec3Dfp(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}


//using fp_float_coords_t = Vec3Dfp;