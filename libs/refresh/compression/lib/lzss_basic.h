#pragma once

#include <cstdint>
#include <vector>
#include <type_traits>
#include <limits>
#include <cassert>

#include <refresh/range_coder/lib/rc.h>

namespace refresh {
	// **********************************************************************************
	// Basic class for LZSS encoding and decosing by range coder
	//template<typename T = uint8_t>
	class lzss_rc_basic
	{
	protected:
		std::vector<uint8_t>& vec;
		refresh::vector_io_stream vios;

		using rc_flag_t = refresh::rc_simple_fixed<refresh::vector_io_stream, 3, 1 << 12, 2>;
		using rc_literal_t = refresh::rc_simple<refresh::vector_io_stream, 1 << 15, 8>;
		using rc_len_t = refresh::rc_simple<refresh::vector_io_stream, 1 << 15, 16>;
		using rc_offset_prefix_t = refresh::rc_simple<refresh::vector_io_stream, 1 << 15, 8>;
		using rc_offset_suffix_t = refresh::rc_simple<refresh::vector_io_stream, 1 << 15, 8>;

		refresh::rc_context_vec_emb<rc_flag_t> dict_flag;
		refresh::rc_context_vec_emb<rc_literal_t> dict_literal;
		refresh::rc_context_vec_emb<rc_len_t> dict_len;
		refresh::rc_context_vec_emb<rc_offset_prefix_t> dict_offset_prefix;
		refresh::rc_context_vec_emb<rc_offset_suffix_t> dict_offset_suffix;

		rc_flag_t* tpl_flag = nullptr;
		rc_literal_t* tpl_literal = nullptr;
		rc_len_t* tpl_len = nullptr;
		rc_offset_prefix_t* tpl_offset_prefix = nullptr;
		rc_offset_suffix_t* tpl_offset_suffix = nullptr;

		refresh::rc_context_t ctx_flag{};
		refresh::rc_context_t ctx_flag_mask = 0xf;

		uint32_t min_len{};
		uint32_t max_len{};
		uint32_t max_offset{};
		uint32_t max_literal_value{};
		uint32_t len_bits{};
		uint32_t offset_bits{};
		uint32_t offset_prefix_bits{};
		uint32_t offset_suffix_bits{};
		uint32_t offset_suffix_mask{};

		template<typename RC>
		void init_tpls(RC& rc, bool compressing)
		{
			tpl_flag = new rc_flag_t(&rc, nullptr, compressing);
			tpl_literal = new rc_literal_t(&rc, max_literal_value + 1, nullptr, compressing);
			tpl_len = new rc_len_t(&rc, max_len - min_len + 1, nullptr, compressing);
			if (offset_prefix_bits)
				tpl_offset_prefix = new rc_offset_prefix_t(&rc, 1 << offset_prefix_bits, nullptr, compressing);
			tpl_offset_suffix = new rc_offset_suffix_t(&rc, 1 << offset_suffix_bits, nullptr, compressing);
		}

	public:
		lzss_rc_basic(std::vector<uint8_t>& vec) :
			vec(vec), vios(vec)
		{
		}

		~lzss_rc_basic()
		{
			if (tpl_flag)
				delete tpl_flag;
			if (tpl_literal)
				delete tpl_literal;
			if (tpl_len)
				delete tpl_len;
			if (tpl_offset_prefix)
				delete tpl_offset_prefix;
			if (tpl_offset_suffix)
				delete tpl_offset_suffix;
		}
	};

	// **********************************************************************************
	//
	// **********************************************************************************
	template<typename T = uint8_t, bool HUGE_DATA = false>
	class lzss_basic
	{
	protected:
		using pos_t = typename std::conditional<HUGE_DATA, uint64_t, uint32_t>::type;

		struct matching_params_t
		{
			uint32_t max_offset;
			uint32_t min_length;
			uint32_t max_length;
			uint32_t max_insert_tries;
			uint32_t max_query_tries;
			bool try_next_pos;

			matching_params_t(uint32_t max_offset = 1 << 15, uint32_t min_length = 3, uint32_t max_length = 1 << 8, uint32_t max_insert_tries = 256, uint32_t max_query_tries = 256, bool try_next_pos = false) :
				max_offset(max_offset),
				min_length(min_length),
				max_length(max_length),
				max_insert_tries(max_insert_tries),
				max_query_tries(max_query_tries),
				try_next_pos(try_next_pos)
			{
			}
		};

		static inline matching_params_t matching_params[] = {
			matching_params_t(1 << 12, 4, 1 << 5, 32, 32, false),				// 0	
			matching_params_t(1 << 13, 4, 1 << 6, 48, 48, false),				// 1
			matching_params_t(1 << 14, 4, 1 << 7, 64, 64, false),				// 2
			matching_params_t(1 << 15, 4, 1 << 8, 64, 64, false),				// 3
			matching_params_t(1 << 16, 4, 1 << 8, 96, 96, false),				// 4
			matching_params_t(1 << 17, 4, 1 << 8, 96, 96, false),				// 5
			matching_params_t(1 << 18, 4, 1 << 8, 128, 128, false),				// 6
			matching_params_t(1 << 19, 4, 1 << 8, 128, 128, true),				// 7
			matching_params_t(1 << 19, 4, 1 << 8, 256, 256, true),				// 8
			matching_params_t(1 << 20, 4, 1 << 8, 512, 512, true),				// 9
			matching_params_t(1 << 21, 4, 1 << 9, 768, 768, true)				// 10
		};

	public:
		lzss_basic()
		{
		}

		static uint32_t get_levels() 
		{
			return (uint32_t)sizeof(matching_params) / sizeof(matching_params[0]);
		}
	};

	// **********************************************************************************
	inline uint32_t lzss_rc_get_levels()
	{
		return lzss_basic<uint8_t, false>::get_levels();
	}
}
