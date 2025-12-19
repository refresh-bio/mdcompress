#pragma once

#include "lzss_basic.h"

namespace refresh {
	// **********************************************************************************
	// Class for LZSS encoding from vector (for testing purposes)
	// **********************************************************************************
	template<typename T = uint32_t>
	class lzss_vector_input
	{
		std::vector<T>& input;
		std::vector<T>::iterator iter;

		const T literal_flag = T(1) << (sizeof(T) * 8 - 1);

		uint32_t min_len{};
		uint32_t max_len{};
		uint32_t max_offset{};
		uint32_t max_literal_value{};
		uint32_t len_bits{};

		T len_mask{};

	public:
		lzss_vector_input(std::vector<T>& input) :
			input(input)
		{
		}

		template<typename U>
		bool get_params(uint32_t& level, U& max_literal_value)
		{
			iter = input.begin();

			if (input.size() < 2)
				return false;

			level = *iter++;
			max_literal_value = *iter++;

			return true;
		}

		void start(uint32_t min_len, uint32_t max_len, uint32_t max_offset)
		{
			this->min_len = min_len;
			this->max_len = max_len;
			this->max_offset = max_offset;

			len_bits = std::bit_width(max_len - min_len);
			len_mask = (T(1) << len_bits) - 1;
		}

		bool is_eod()
		{
			return iter == input.end();
		}

		bool is_literal()
		{
			return (*iter & literal_flag) != 0;
		}

		// Encode a single literal byte
		template<typename U>
		void decode_literal(U& literal)
		{
			literal = U((*iter++ & ~literal_flag));
		}

		// Encode a match with given offset and length
		void decode_match(uint32_t& offset, uint32_t& length)
		{
			length = (*iter & len_mask) + min_len;
			offset = *iter++ >> len_bits;
		}

		void complete()
		{
		}
	};

	// **********************************************************************************
	//
	// **********************************************************************************
	class lzss_rc_input : public lzss_rc_basic
	{
		refresh::rc_decoder<refresh::vector_io_stream> rcd;

		bool eod_noticed = false;

	public:
		lzss_rc_input(std::vector<uint8_t>& vec) :
			lzss_rc_basic(vec),
			rcd(lzss_rc_basic::vios)
		{
		}

		bool get_params(uint32_t& level, uint32_t& max_literal_value)
		{
			vios.restart_read();

			rcd.start();

			level = rcd.get_cumulative_freq(256);
			rcd.update_frequency(1, level, 256);
			max_literal_value = 0;

			for (uint32_t i = 0; i < 4; ++i)
			{
				uint32_t x = rcd.get_cumulative_freq(256);
				rcd.update_frequency(1, x, 256);
				max_literal_value |= (x << (i * 8));
			}

			this->max_literal_value = max_literal_value;
			eod_noticed = false;

			ctx_flag = 0;

			return true;
		}

		void start(uint32_t min_len, uint32_t max_len, uint32_t max_offset)
		{
			this->min_len = min_len;
			this->max_len = max_len;
			this->max_offset = max_offset;

			len_bits = std::bit_width(max_len - min_len);
			offset_bits = std::bit_width(max_offset - 1);

			if (offset_bits <= 8)
			{
				offset_prefix_bits = 0;
				offset_suffix_bits = offset_bits;
			}
			else
			{
				offset_prefix_bits = 8;
				offset_suffix_bits = offset_bits - offset_prefix_bits;
			}

			offset_suffix_mask = (1 << offset_suffix_bits) - 1;

			init_tpls(rcd, false);
		}

		bool is_eod()
		{
			return eod_noticed;
		}

		bool is_literal()
		{
			auto flag = rc_find_context(dict_flag, ctx_flag, *tpl_flag)->decode();

			if (flag == 2)
			{
				eod_noticed = true;
				return false;
			}

			ctx_flag = ((ctx_flag << 1) | (uint8_t)flag) & ctx_flag_mask;
			return flag == 0;
		}

		// Decode a single literal
		template<typename U>
		void decode_literal(U& literal)
		{
			literal = rc_find_context(dict_literal, 0, *tpl_literal)->decode();
		}

		// Encode a match with given offset and length
		void decode_match(uint32_t& offset, uint32_t& length)
		{
			if (offset_prefix_bits)
			{
				uint32_t offset_prefix = rc_find_context(dict_offset_prefix, 0, *tpl_offset_prefix)->decode();
				uint32_t offset_suffix;

				if (offset_prefix <= 2)
					offset_suffix = rc_find_context(dict_offset_suffix, offset_prefix, *tpl_offset_suffix)->decode();
				else
				{
					offset_suffix = rcd.get_cumulative_freq(offset_suffix_mask + 1);
					rcd.update_frequency(1, offset_suffix, offset_suffix_mask + 1);
				}

				offset = (offset_prefix << offset_suffix_bits) | offset_suffix;
			}
			else
			{
				offset = rc_find_context(dict_offset_suffix, 0, *tpl_offset_suffix)->decode();
			}

			offset += 1;		// offset is stored as (offset-1)

			length = rc_find_context(dict_len, 0, *tpl_len)->decode() + min_len;
		}

		void complete()
		{
			rcd.complete();
		}
	};

	// **********************************************************************************
	//
	// **********************************************************************************
	template<typename T = uint8_t>
	class lzss_decompressor : public lzss_basic<T, false>
	{
		using pos_t = typename lzss_basic<T, false>::pos_t;

	public:
		// **********************************************************************************
		lzss_decompressor() : lzss_basic<T, false>()
		{
		}

		// **********************************************************************************
		template<typename EC_TYPE>
		bool decompress(std::vector<T>& raw, EC_TYPE& decoder)
		{
			uint32_t level;
			uint32_t max_literal_value;

			raw.clear();

			if (!decoder.get_params(level, max_literal_value))
				return false;

			uint32_t min_length = lzss_basic<T, false>::matching_params[level].min_length;
			uint32_t max_length = lzss_basic<T, false>::matching_params[level].max_length;
			uint32_t max_offset = lzss_basic<T, false>::matching_params[level].max_offset;

			decoder.start(min_length, max_length, max_offset);

			while (!decoder.is_eod())
			{
				if (decoder.is_literal())
				{
					T literal;
					decoder.decode_literal(literal);

					if (literal > max_literal_value)
						return false;

					raw.push_back(literal);
				}
				else if (decoder.is_eod())
					break;
				else
				{
					uint32_t offset, length;

					decoder.decode_match(offset, length);

					if (offset == 0 || offset > max_offset || length < min_length || length > max_length)
						return false;

					pos_t start_pos = (pos_t)raw.size() - offset;

					for (uint32_t i = 0; i < length; ++i)
						raw.push_back(raw[start_pos + i]);
				}
			}

			return true;
		}
	};

	// **********************************************************************************
	inline bool lzss_rc_decompress(std::vector<uint8_t>& packed, std::vector<uint8_t>& raw)
	{
		refresh::lzss_rc_input ed(packed);
		refresh::lzss_decompressor lzss;

		return lzss.decompress(raw, ed);
	}
}
