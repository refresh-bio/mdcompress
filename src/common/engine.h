#pragma once

#include "utils.h"
#include "frame.h"
#include "../libs/refresh/range_coder/lib/rc.h"
#include "vec3d.h"
#include "../libs/refresh/fp_float/lib/fp_float.h"
#include <bit>
#include "../app/params.h"

// *******************************************************************************
template<typename T>
class CEngineState
{
	friend class CEngineAbstract;
	friend class CEngineCompress;
	friend class CEngineDecompress;

	// ***** Range Coder
	//    static const uint32_t prefix_max_bits = 13;                     // !!! max: 15        
	static const uint32_t prefix_max_bits = 14;                     // !!! max: 15        
	static const uint32_t prefix_short_max_bits = 6;
//	static const uint32_t prefix_short_max_bits = 9;
	//    static const uint32_t prefix_model_size = prefix_max_bits + 2;
	static const uint32_t prefix_model_size = prefix_max_bits + 2;

	using model_binary_t = refresh::rc_simple_fixed<refresh::vector_io_stream, 2, 1 << 10, 1>;
	using model_prefix_t = refresh::rc_simple_fixed<refresh::vector_io_stream, prefix_model_size, 1 << 15, 2>;
	using model_suffix_t = refresh::rc_simple<refresh::vector_io_stream, 1 << 14, 1>;

	refresh::rc_context_vec_emb<model_binary_t> dict_binary;
	refresh::rc_context_vec_emb<model_prefix_t> dict_cart_prefix;
	refresh::rc_context_vec<model_suffix_t> dict_cart_suffix;

	model_binary_t* tpl_binary = nullptr;

	model_prefix_t* tpl_prefix = nullptr;
	std::vector<model_suffix_t> tpl_suffix;

	std::vector<uint8_t> v_rc;
	refresh::vector_io_stream* vios = nullptr;
	T* rc = nullptr;

	void clear_rc_dict()
	{
		dict_binary.clear();
		dict_cart_prefix.clear();
		dict_cart_suffix.clear();

		v_rc.clear();
		vios->restart_read();
	}

	template<typename S> 
	void delete_var(S*& var)
	{
		if (var)
		{
			delete var;
			var = nullptr;
		}
	}

	void init_rc_types()
	{
		vios = new refresh::vector_io_stream(v_rc);

		rc = new T(*vios);
	}

	void prepare_rc_tpl(bool is_compressing)
	{
		if (!tpl_binary)
			tpl_binary = new model_binary_t(rc, nullptr, is_compressing);

		if (!tpl_prefix)
			tpl_prefix = new model_prefix_t(rc, nullptr, is_compressing);

		if (tpl_suffix.empty())
		{
			for (uint32_t i = 0; i <= prefix_short_max_bits; ++i)
				tpl_suffix.push_back(model_suffix_t(rc, 1 << i, nullptr, is_compressing));
			tpl_suffix.push_back(model_suffix_t(rc, 256, nullptr, is_compressing));
		}
	}

	void release_rc_types()
	{
		delete_var(vios);
		delete_var(rc);

		delete_var(tpl_binary);
		delete_var(tpl_prefix);
		tpl_suffix.clear();
	}

public:
	CEngineState(bool is_compressing)
	{
		init_rc_types();
		prepare_rc_tpl(is_compressing);
	}

	~CEngineState() 
	{
		release_rc_types();
	}
};

// *******************************************************************************
class CEngineAbstract
{
protected:
	using my_float = Vec3Dfp::my_float;

	compression_features_t compression_features;

	template<typename T> uint32_t no_bits(T x)
	{
#ifdef __cpp_lib_int_pow2
		return std::bit_width(x);
#else
		uint32_t r;

		for (r = 0; x; ++r)
			x >>= 1;

		return r;
#endif
	}

	template<typename T> int64_t ipow2(T a)
	{
		return (int64_t)a * (int64_t)a;
	}

	// Function sqrt()
/*	template <typename T>
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
*/	
	// Predictors
	// *******************************************************************************
	std::pair<Vec3Dfp, Vec3Dfp> get_q_both(const Vec3Dfp& b1, const Vec3Dfp& b2, const Vec3Dfp& b3, const my_float& dist_q1_sq, const my_float& dist_q2_sq, const my_float& dist_q3_sq, bool ref_left_handed)
	{
		Vec3Dfp b21 = b2 - b1;
		Vec3Dfp b31 = b3 - b1;

		my_float d = vector_length(b21);

		Vec3Dfp ex = b21 / d;

		my_float i = dot_product(b31, ex);

		Vec3Dfp temp = b31 - (i * ex);

		my_float j_sq = vector_length2(temp);
		my_float j = sqrt(j_sq);
		Vec3Dfp ey = temp / j;

		Vec3Dfp ez = cross_product(ex, ey);

		my_float d_sq = pow2(d);
		my_float i_sq = pow2(i);
//		my_float j_sq = pow2(j);

		my_float x = ((dist_q1_sq - dist_q2_sq + d_sq) / d) >> 1;

		my_float numerator_y = dist_q1_sq - dist_q3_sq + i_sq + j_sq - (i * x * 2);
		my_float y = (numerator_y / j) >> 1;

		my_float z2 = dist_q1_sq - pow2(x) - pow2(y);

		Vec3Dfp q_base = b1 + (x * ex) + (y * ey);

		if (z2.is_positive())
		{
			my_float z = sqrt(z2);
			ez = z * ez;

			auto r1 = q_base + ez;
			auto r2 = q_base - ez;

			auto r1_ok = !r1.is_wrong();
			auto r2_ok = !r2.is_wrong();

			if (r1_ok && r2_ok)
			{
				if (compression_features.molecule_orientation == compression_features_t::molecule_orientation_t::predicted)
				{
					if (is_left_handed(b1, b2, b3, r1) == ref_left_handed)
						return std::make_pair(r2, r1);
					else
						return std::make_pair(r1, r2);
				}
				else
					return std::make_pair(r1, r2);
			}

			if (r1_ok)
				return std::make_pair(r1, r1);

			if (r2_ok)
				return std::make_pair(r2, r2);

			return std::make_pair(q_base, q_base);
		}
		else
			return std::make_pair(q_base, q_base);
	}

	// *******************************************************************************
	Vec3Dfp get_q_new2(const Vec3Dfp& a0, const Vec3Dfp& b0, const Vec3Dfp& p0, const Vec3Dfp& a1, const Vec3Dfp& b1, const Vec3Dfp& p1, const Vec3Dfp& q0)
	{
		Vec3Dfp e1 = normalize_vector(a0 - p0);
		Vec3Dfp vb0 = b0 - p0;
		Vec3Dfp::my_float dot0 = dot_product(vb0, e1);
		Vec3Dfp proj0 = dot0 * e1;
		Vec3Dfp e2 = normalize_vector(vb0 - proj0);
		Vec3Dfp e3 = cross_product(e1, e2);

		Vec3Dfp f1 = normalize_vector(a1 - p1);
		Vec3Dfp vb1 = b1 - p1;
		Vec3Dfp::my_float dot1 = dot_product(vb1, f1);
		Vec3Dfp proj1 = dot1 * f1;
		Vec3Dfp f2 = normalize_vector(vb1 - proj1);
		Vec3Dfp f3 = cross_product(f1, f2);

		auto R_00 = f1.x * e1.x + f2.x * e2.x + f3.x * e3.x;
		auto R_01 = f1.x * e1.y + f2.x * e2.y + f3.x * e3.y;
		auto R_02 = f1.x * e1.z + f2.x * e2.z + f3.x * e3.z;

		auto R_10 = f1.y * e1.x + f2.y * e2.x + f3.y * e3.x;
		auto R_11 = f1.y * e1.y + f2.y * e2.y + f3.y * e3.y;
		auto R_12 = f1.y * e1.z + f2.y * e2.z + f3.y * e3.z;

		auto R_20 = f1.z * e1.x + f2.z * e2.x + f3.z * e3.x; 
		auto R_21 = f1.z * e1.y + f2.z * e2.y + f3.z * e3.y;
		auto R_22 = f1.z * e1.z + f2.z * e2.z + f3.z * e3.z;

		auto dx = q0.x - p0.x;
		auto dy = q0.y - p0.y;
		auto dz = q0.z - p0.z;

		return Vec3Dfp(
			p1.x + R_00 * dx + R_01 * dy + R_02 * dz, 
			p1.y + R_10 * dx + R_11 * dy + R_12 * dz, 
			p1.z + R_20 * dx + R_21 * dy + R_22 * dz
		);
	}

	// *******************************************************************************
	std::pair<int, int> pred_z1(const int64_t d2_OH, const int x0, const int y0, const int z0, const int x1, const int y1)
	{
		int64_t delta = d2_OH - ipow2(x0 - x1) - ipow2(y0 - y1);

		if (delta < 0)
			return std::make_pair(z0, z0);

		int is_delta = static_cast<int64_t>(isqrt(static_cast<uint64_t>(delta)));

		return std::make_pair(z0 - is_delta, z0 + is_delta);
	}

	// *******************************************************************************
	int pred_y2_old(int64_t d2_OH, int64_t d2_HH, int x0, int y0, int z0, int x1, int y1, int z1, int x2)
	{
		double C = (double)(ipow2(x0 - x2) - d2_OH);
		double D = (double)(ipow2(x1 - x2) - d2_HH);
		double E = (double)(ipow2(y0) - ipow2(y1) + ipow2(z0) - ipow2(z1) + C - D);

		if (z1 - z0 == 0)
		{
			if (y1 - y0 == 0)
				return y0;

			return iround<int>(-ipow2(y0) + ipow2(y1) - ipow2(z0) + ipow2(z1) - C + D) / (2 * y1 - 2 * y0);
		}

		double G = double(y0 - y1) / (z1 - z0);
		double F = E / (2 * z1 - 2 * z0);
		double K = 1 + G * G;
		double L = -2 * y1 - 2 * G * F - 2 * z1 * G;
		double M = ipow2(y1) + ipow2(z1) + 2 * z1 * F + F * F + D;
		double delta = L * L - 4 * K * M;

		if (delta < 0)
			return iround<int>(-L / (2 * K));

		double s_delta = sqrt(delta);
		int iy2_a = iround<int>((-L + s_delta) / (2 * K));
		int iy2_b = iround<int>((-L - s_delta) / (2 * K));

		int64_t pred_z2a = pred_z2(d2_OH, d2_HH, x0, y0, z0, x1, y1, z1, x2, iy2_a);
		int64_t pred_z2b = pred_z2(d2_OH, d2_HH, x0, y0, z0, x1, y1, z1, x2, iy2_b);

		int64_t dx2_02 = ipow2(x2 - x0);
		int64_t d2_OH_a = dx2_02 + ipow2(iy2_a - y0) + ipow2(pred_z2a - z0);
		int64_t d2_OH_b = dx2_02 + ipow2(iy2_b - y0) + ipow2(pred_z2b - z0);

		// Choose the candidate with the smaller error
		if (std::abs(d2_OH - d2_OH_a) < std::abs(d2_OH - d2_OH_b))
			return iy2_a;
		else
			return iy2_b;
	}

	// *******************************************************************************
	int pred_y2(const int64_t d2_OH, const int64_t d2_HH, const int x0, const int y0, const int z0, const int x1, const int y1, const int z1, const int x2)
	{
		const int64_t C = ipow2(x0 - x2) - d2_OH;
		const int64_t D = ipow2(x1 - x2) - d2_HH;
		const int64_t E = ipow2(y0) - ipow2(y1) + ipow2(z0) - ipow2(z1) + C - D;

		const int64_t dz10 = z1 - z0;
		const int64_t dy01 = y0 - y1;

		if (dz10 == 0)
		{
			if (dy01 == 0)
				return y0;				// Very rare case, so we do not care about it

			return idiv<int64_t>(E, 2 * dy01);
		}

		my_float G;

		bool is_G_inversed = abs(dy01) < abs(dz10);
		if (is_G_inversed)
			G = my_float(dz10) / my_float(dy01);
		else
			G = my_float(dy01) / my_float(dz10);

		my_float F = my_float(E) / my_float(2 * dz10);
		int64_t extension = 16;
		my_float K_extended = 1 << extension;
		my_float pow2_G = pow2(G);
		if(is_G_inversed)
			K_extended += my_float(1 << extension) / pow2_G;
		else
			K_extended += pow2_G << extension;
		my_float L = -y1;
		
		if(is_G_inversed)
			L -= (F / G) + (my_float(z1) / G);
		else
			L -= (G * F) + (G * my_float(z1));
		L <<= 1;

		my_float M = my_float(ipow2(y1) + ipow2(z1)) + (F * (2 * z1)) + pow2(F) + D;
		my_float M_4 = M << 2;
		my_float delta = pow2(L) - M_4;
		if (is_G_inversed)
			delta -= M_4 / pow2_G;
		else
			delta -= M_4 * pow2_G;

		my_float KK = K_extended << 1;

		if (delta.is_negative())
		{
			my_float y2 = ((-L) << extension) / KK;

			if (y2.is_wrong())
				return y0;
			return (int) int_round_or_0(y2);
		}

		my_float s_delta = sqrt(delta);
//		my_float y2_a = round(((-L + s_delta) << extension) / KK);
//		my_float y2_b = round(((-L - s_delta) << extension) / KK);
		my_float y2_a = ((-L + s_delta) << extension) / KK;
		my_float y2_b = ((-L - s_delta) << extension) / KK;

		if (y2_a.is_wrong() || y2_b.is_wrong())
			return y0;

//		int iy2_a = (int) y2_a;
//		int iy2_b = (int) y2_b;
		int iy2_a = (int) int_round_or_0(y2_a);
		int iy2_b = (int)int_round_or_0(y2_b);

		int64_t pred_z2a = pred_z2(d2_OH, d2_HH, x0, y0, z0, x1, y1, z1, x2, iy2_a);
		int64_t pred_z2b = pred_z2(d2_OH, d2_HH, x0, y0, z0, x1, y1, z1, x2, iy2_b);

		int64_t dx2_02 = ipow2(x2 - x0);
		int64_t d2_OH_a = dx2_02 + ipow2(iy2_a - y0) + ipow2(pred_z2a - z0);
		int64_t d2_OH_b = dx2_02 + ipow2(iy2_b - y0) + ipow2(pred_z2b - z0);

		// Choose the candidate with the smaller error
		if (std::abs(d2_OH - d2_OH_a) < std::abs(d2_OH - d2_OH_b))
			return iy2_a;
		else
			return iy2_b;
	}

	// *******************************************************************************
	int pred_x2(const int64_t d_OH, const int64_t d_HH, const int64_t d2_OH, const int64_t d2_HH, const int x0, const int y0, const int z0, const int x1, const int y1, const int z1)
	{
		Vec3Dfp p0(x0, y0, z0);
		Vec3Dfp p1(x1, y1, z1);

		auto d_vec = p1 - p0;
		my_float d = vector_length(d_vec);

		if (d > d_OH + d_HH)
			return x0; 

		if (d < iabs(d_OH - d_HH))
			return x0;

		my_float a = (pow2(d) + d2_OH - d2_HH) / (d * 2);

		auto d_vec_normalized = normalize_vector(d_vec);
		auto out_center = p0 + a * d_vec_normalized;

		if (out_center.is_wrong())
			return x0; 

		return iround<int>(out_center.x);
	}

	// *******************************************************************************
	std::pair<int, int> pred_y2_new(const int64_t d2_OH, const int64_t d2_HH, const int x0, const int y0, const int z0, const int x1, const int y1, const int z1, const int x2)
	{
		const int64_t C = ipow2(x0 - x2) - d2_OH;
		const int64_t D = ipow2(x1 - x2) - d2_HH;
		const int64_t E = ipow2(y0) - ipow2(y1) + ipow2(z0) - ipow2(z1) + C - D;

		const int64_t dz10 = z1 - z0;
		const int64_t dy01 = y0 - y1;

		if (dz10 == 0)
		{
			if (dy01 == 0)
				return std::make_pair(y0, y0);				// Very rare case, so we do not care about it

			auto a = idiv<int64_t>(E, 2 * dy01);
			return std::make_pair((int)a, (int)a);
		}

		my_float G;

		bool is_G_inversed = abs(dy01) < abs(dz10);
		if (is_G_inversed)
			G = my_float(dz10) / my_float(dy01);
		else
			G = my_float(dy01) / my_float(dz10);

		my_float F = my_float(E) / my_float(2 * dz10);
		int64_t extension = 16;
		my_float K_extended = 1 << extension;
		my_float pow2_G = pow2(G);
		if(is_G_inversed)
			K_extended += my_float(1 << extension) / pow2_G;
		else
			K_extended += pow2_G << extension;
		my_float L = -y1;
		
		if(is_G_inversed)
			L -= (F / G) + (my_float(z1) / G);
		else
			L -= (G * F) + (G * my_float(z1));
		L <<= 1;

		my_float M = my_float(ipow2(y1) + ipow2(z1)) + (F * (2 * z1)) + pow2(F) + D;
		my_float M_4 = M << 2;
		my_float delta = pow2(L) - M_4;
		if (is_G_inversed)
			delta -= M_4 / pow2_G;
		else
			delta -= M_4 * pow2_G;

		my_float KK = K_extended << 1;

		if (delta.is_negative())
		{
			my_float y2 = ((-L) << extension) / KK;

			if (y2.is_wrong())
				return std::make_pair(y0, y0);
			return std::make_pair((int)int_round_or_0(y2), (int)int_round_or_0(y2));
		}

		my_float s_delta = sqrt(delta);
//		my_float y2_a = round(((-L + s_delta) << extension) / KK);
//		my_float y2_b = round(((-L - s_delta) << extension) / KK);
		my_float y2_a = ((-L + s_delta) << extension) / KK;
		my_float y2_b = ((-L - s_delta) << extension) / KK;

		if (y2_a.is_wrong() || y2_b.is_wrong())
			return std::make_pair(y0, y0);

		return std::make_pair((int)int_round_or_0(y2_a), (int)int_round_or_0(y2_b));
	}
	
	// *******************************************************************************
	int pred_z2(const int64_t d2_OH, const int64_t d2_HH, const int x0, const int y0, const int z0, const int x1, const int y1, const int z1, const int x2, const int y2)
	{
		const int64_t delta = d2_OH - ipow2(x0 - x2) - ipow2(y0 - y2);

		if (delta < 0)
			return z0;

		int is_delta = static_cast<int>(isqrt(static_cast<uint64_t>(delta)));

		int iz_a = z0 - is_delta;
		int iz_b = z0 + is_delta;

		int64_t d2_common = ipow2(x2 - x1) + ipow2(y2 - y1);
		int64_t d2_first = d2_common + ipow2(iz_a - z1);
		int64_t d2_second = d2_common + ipow2(iz_b - z1);

		return (abs(d2_HH - d2_first) < abs(d2_HH - d2_second)) ? iz_a : iz_b;
	}

	// *******************************************************************************
	icoord_t history_predictor(icoord_t& a)
	{
		return a;
	}

	// *******************************************************************************
	icoord_t history_predictor(icoord_t& a, icoord_t& b)
	{
		return icoord_t(2 * a.x - b.x, 2 * a.y - b.y, 2 * a.z - b.z);
	}

	// *******************************************************************************
	icoord_t history_predictor(icoord_t& a, icoord_t& b, icoord_t& c)
	{
		return icoord_t(3 * (a.x - b.x) + c.x, 3 * (a.y - b.y) + c.y, 3 * (a.z - b.z) + c.z);
	}

public:
	CEngineAbstract()
	{
	}
	
	virtual ~CEngineAbstract() = default;

	void set_compression_features(const compression_features_t& cf)
	{
		compression_features = cf;
	}

	// *******************************************************************************
	bool is_left_handed(const Vec3Dfp& p0, const Vec3Dfp& p1, const Vec3Dfp& p2, const Vec3Dfp& q) const
	{
		const auto v1 = p1 - p0;
		const auto v2 = p2 - p0;
		const auto vq = q - p0;

		auto normal_vector = cross_product(v1, v2);
		auto mixed_product = dot_product(normal_vector, vq);

		return mixed_product.is_negative();
	}
	// *******************************************************************************
	bool is_left_handed(const icoord_t& p0, const icoord_t& p1, const icoord_t& p2, const icoord_t& q) const
	{
		return is_left_handed(Vec3Dfp(p0.x, p0.y, p0.z), Vec3Dfp(p1.x, p1.y, p1.z), Vec3Dfp(p2.x, p2.y, p2.z), Vec3Dfp(q.x, q.y, q.z));
	}
};
