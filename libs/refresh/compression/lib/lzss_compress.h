#pragma once

#include <cstdio>
#include <iostream>
#include <bit>

#include <refresh/string_operations/lib/string_operations.h>

#include "lzss_basic.h"

namespace refresh {
	// **********************************************************************************
	// Class for LZSS encoding to stdout (for testing purposes)
	// **********************************************************************************
	class lzss_ec_stdout
	{
	public:
		lzss_ec_stdout() = default;

		void start(size_t level, size_t min_len, size_t max_len, size_t max_offset, size_t max_literal_value)
		{
		}

		// Encode a single literal byte
		void encode_literal(uint8_t literal)
		{
			std::cout << literal << " ";
		}

		// Encode a match with given offset and length
		void encode_match(unsigned offset, unsigned length)
		{
			// For demonstration, we print the match as (offset,length)
			std::cout << "(" << offset << "," << length << ") ";
		}

		void complete()
		{
		}
	};

	// **********************************************************************************
	// Class for LZSS encoding to vector (for testing purposes)
	// **********************************************************************************
	template<typename T = uint32_t>
	class lzss_vector_output
	{
		std::vector<T>& output;
		const T literal_flag = T(1) << (sizeof(T) * 8 - 1);

		uint32_t min_len{};
		uint32_t max_len{};
		uint32_t max_offset{};
		uint32_t max_literal_value{};
		uint32_t len_bits{};

	public:
		lzss_vector_output(std::vector<T>& output) :
			output(output)
		{
		}

		void start(uint32_t level, uint32_t min_len, uint32_t max_len, uint32_t max_offset, T max_literal_value)
		{
			this->min_len = min_len;
			this->max_len = max_len;
			this->max_offset = max_offset;
			this->max_literal_value = max_literal_value;

			len_bits = std::bit_width(max_len - min_len);

			output.push_back(T(level));
			output.push_back(T(max_literal_value));
		}

		// Encode a single literal byte
		template<typename U>
		void encode_literal(U literal)
		{
			output.push_back(T(literal) + literal_flag);
		}

		// Encode a match with given offset and length
		void encode_match(uint32_t offset, uint32_t length)
		{
			output.push_back(T(offset << len_bits) | T(length - min_len));
		}

		void complete()
		{
		}
	};

	// **********************************************************************************
	// Class for LZSS encoding by range coder
	// **********************************************************************************
	class lzss_rc_output : public lzss_rc_basic
	{
		refresh::rc_encoder<refresh::vector_io_stream> rce;

		using rc_flag_t = typename lzss_rc_basic::rc_flag_t;
		using rc_literal_t = typename lzss_rc_basic::rc_literal_t;
		using rc_len_t = typename lzss_rc_basic::rc_len_t;
		using rc_offset_prefix_t = typename lzss_rc_basic::rc_offset_prefix_t;
		using rc_offset_suffix_t = typename lzss_rc_basic::rc_offset_suffix_t;

	public:
		lzss_rc_output(std::vector<uint8_t>& vec) :
			lzss_rc_basic(vec),
			rce(lzss_rc_basic::vios)
		{
			vios.clear();
		}

		~lzss_rc_output() = default;

		void start(uint32_t level, uint32_t min_len, uint32_t max_len, uint32_t max_offset, uint32_t max_literal_value)
		{
			this->min_len = min_len;
			this->max_len = max_len;
			this->max_offset = max_offset;
			this->max_literal_value = max_literal_value;

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

			init_tpls(rce, true);

			rce.start();

			rce.encode_frequency(1, level, 256);

			uint32_t x = max_literal_value;
			for (uint32_t i = 0; i < 4; ++i)
			{
				rce.encode_frequency(1, x & 0xff, 256);
				x >>= 8;
			}

			ctx_flag = 0;
		}

		// Encode a single literal byte
		template<typename U>
		void encode_literal(U literal)
		{
			rc_find_context(dict_flag, ctx_flag, *tpl_flag)->encode(0);
			ctx_flag = ((ctx_flag << 1) | 0) & ctx_flag_mask;

			rc_find_context(dict_literal, 0, *tpl_literal)->encode(literal);
		}

		// Encode a match with given offset and length
		void encode_match(uint32_t offset, uint32_t length)
		{
			rc_find_context(dict_flag, ctx_flag, *tpl_flag)->encode(1);
			ctx_flag = ((ctx_flag << 1) | 1) & ctx_flag_mask;

			offset -= 1;		// offset is stored as (offset-1)

			if (offset_prefix_bits)
			{
				uint32_t offset_prefix = offset >> offset_suffix_bits;
				uint32_t offset_suffix = offset & offset_suffix_mask;

				rc_find_context(dict_offset_prefix, 0, *tpl_offset_prefix)->encode(offset_prefix);

				if (offset_prefix <= 2)
					rc_find_context(dict_offset_suffix, offset_prefix, *tpl_offset_suffix)->encode(offset_suffix);
				else
					rce.encode_frequency(1, offset_suffix, offset_suffix_mask + 1);
			}
			else
			{
				rc_find_context(dict_offset_suffix, 0, *tpl_offset_suffix)->encode(offset);
			}

			rc_find_context(dict_len, 0, *tpl_len)->encode(length - min_len);
		}

		void complete()
		{
			rc_find_context(dict_flag, ctx_flag, *tpl_flag)->encode(2);			// EOD
			rce.complete();
		}
	};

	// **********************************************************************************
	//
	// **********************************************************************************
	template<typename T = uint8_t, bool HUGE_DATA = false>
	class lzss_compressor : public lzss_basic<T, HUGE_DATA>
	{
		using pos_t = typename lzss_basic<T, HUGE_DATA>::pos_t;

		const pos_t empty_entry = std::numeric_limits<pos_t>::max();
		const pos_t removed_entry = empty_entry - 1;
		const double ht_load_factor = 0.5;
		size_t ht_size{};
		size_t ht_mask{};

		std::vector<pos_t> ht, ht_aux;

		// Clear and (re)allocate the hash table
		void clear_ht(size_t max_offset)
		{
			ht_size = std::bit_ceil(size_t(max_offset * 2 / ht_load_factor));
			ht_mask = ht_size - 1;

			ht.clear();
			ht.resize(ht_size, empty_entry);
		}

		// Refresh the hash table - remove entries that are out of the sliding window
		void refresh_ht(const T* input, size_t pos, size_t max_offset, size_t min_length, uint32_t max_insert_tries)
		{
			ht_aux.clear();
			ht_aux.swap(ht);
			ht.resize(ht_size, empty_entry);

			for (pos_t i = 0; i < (pos_t)ht_aux.size(); ++i)
			{
				auto e = ht[i];
				if (e != empty_entry && e + max_offset >= pos)
				{
					auto h = hasher(input + e, min_length);
					insert_ht(h, e, max_offset, max_insert_tries);
				}
			}
		}

		void insert_ht(size_t h, pos_t pos, uint32_t max_offset, uint32_t max_insert_tries)
		{
			if (h + max_insert_tries <= ht_size)
			{
				uint32_t h_max = h + max_insert_tries;
				for (; h < h_max; ++h)
					if (ht[h] >= removed_entry)
					{
						ht[h] = pos;
						return;
					}

				return;
			}
			else
			{
				uint32_t h_max = (h + max_insert_tries) & ht_mask;

				for (; h < ht_size; ++h)
					if (ht[h] >= removed_entry)
					{
						ht[h] = pos;
						return;
					}

				for (h = 0; h < h_max; ++h)
					if (ht[h] >= removed_entry)
					{
						ht[h] = pos;
						return;
					}
			}
		}

		// A simple but effective hash function for short sequences (up to 5 bytes)
		size_t hasher(const T* p, uint32_t len) const
		{
			size_t h = 0;

			switch (len)
			{
			case 5:	h ^= size_t(p[4]) * size_t(0xff51afd7ed558ccdL);
			case 4:	h ^= size_t(p[3]) * size_t(0xc4ceb9fe1a85ec53L);
			case 3:	h ^= size_t(p[2]) * size_t(0x9e3779b97f4a7c15L);
			case 2:	h ^= size_t(p[1]) * size_t(0x94d049bb133111ebL);
			case 1:	h ^= size_t(p[0]) * size_t(0x2545f4914f6cdd1dL);
				//					h = (h ^ (h >> 33)) * size_t(0xff51afd7ed558ccdL);
				//					h = (h ^ (h >> 33)) * size_t(0xc4ceb9fe1a85ec53L);
				//					h = h ^ (h >> 33);
				break;
			default:
				assert(false);
			}

			return h & ht_mask;
		}

		std::pair<uint32_t, uint32_t> find_match(const std::vector<T>& raw, size_t h, pos_t pos, uint32_t max_offset, uint32_t max_length, uint32_t max_tries, uint32_t& possible_insert_pos)
		{
			uint32_t best_len = 0;
			uint32_t best_pos = 0;

			possible_insert_pos = ~uint32_t(0);

			uint32_t effective_max_length = std::min<uint32_t>(max_length, (pos_t)raw.size() - pos);

			for (uint32_t tries = 0; tries < max_tries; ++tries, h = (h + 1) & ht_mask)
			{
				auto e = ht[h];
				if (e >= removed_entry)
				{
					if (possible_insert_pos == ~uint32_t(0))
						possible_insert_pos = h;

					if (e == empty_entry)
						break;
					else
						continue;
				}

				if (pos - e > max_offset)
				{
					if (possible_insert_pos == ~uint32_t(0))
						possible_insert_pos = h;

					ht[h] = removed_entry;

					continue;
				}

				if (best_len > 1 && raw[e + best_len - 1] != raw[pos + best_len - 1])
					continue;

				uint32_t match_len = (uint32_t)refresh::matching_length(&raw[e], &raw[pos], effective_max_length);

				if (match_len >= best_len)
				{
					if (match_len > best_len)
					{
						best_len = match_len;
						best_pos = e;

						if (best_len == effective_max_length)
							break;
					}
					else if (e > best_pos)
						best_pos = e;
				}
			}

			return { best_pos, best_len };
		}

	public:
		// **********************************************************************************
		lzss_compressor() : lzss_basic<T, HUGE_DATA>()
		{
		}

		// **********************************************************************************
		template<typename EC_TYPE>
		bool compress(const std::vector<T>& raw, EC_TYPE& encoder, uint32_t level, T max_literal_val = std::numeric_limits<T>::max())
		{
			if (level >= lzss_basic<T, HUGE_DATA>::get_levels())
				level = lzss_basic<T, HUGE_DATA>::get_levels() - 1;

			uint32_t min_length = lzss_basic<T, HUGE_DATA>::matching_params[level].min_length;
			uint32_t max_length = lzss_basic<T, HUGE_DATA>::matching_params[level].max_length;
			uint32_t max_offset = lzss_basic<T, HUGE_DATA>::matching_params[level].max_offset;
			uint32_t max_insert_tries = lzss_basic<T, HUGE_DATA>::matching_params[level].max_insert_tries;
			uint32_t max_query_tries = lzss_basic<T, HUGE_DATA>::matching_params[level].max_query_tries;
			bool try_next_pos = lzss_basic<T, HUGE_DATA>::matching_params[level].try_next_pos;

			pos_t next_refresh = max_offset;
			pos_t pos;

			clear_ht(max_offset);

			encoder.start(level, min_length, max_length, max_offset, max_literal_val);

			for (pos = 0; pos + min_length < (pos_t)raw.size();)
			{
				pos_t best_pos = 0;
				pos_t best_len = 0;
//				uint32_t tries = 0;

				auto h = hasher(&raw[pos], min_length);

				uint32_t possible_insert_pos;

				std::tie(best_pos, best_len) = find_match(raw, h, pos, max_offset, max_length, max_query_tries, possible_insert_pos);

				if (try_next_pos && best_len >= min_length && pos + min_length + 1 < raw.size())
				{
					auto h2 = hasher(&raw[pos + 1], min_length);
					uint32_t best_pos2, best_len2, possible_insert_pos2;
					std::tie(best_pos2, best_len2) = find_match(raw, h2, pos + 1, max_offset, max_length, max_query_tries, possible_insert_pos2);

					if (best_len2 > best_len)
					{
						encoder.encode_literal(raw[pos]);

						if (possible_insert_pos != ~uint32_t(0))
							h = possible_insert_pos;

						insert_ht(h, pos++, max_offset, max_insert_tries);

						best_pos = best_pos2;
						best_len = best_len2;
						possible_insert_pos = possible_insert_pos2;
					}
				}

				if (best_len >= min_length)
				{
					encoder.encode_match(pos - best_pos, best_len);

					// Insert all sequences in the matched string to the hash table
					if (possible_insert_pos != ~uint32_t(0))
						h = possible_insert_pos;

					insert_ht(h, pos, max_offset, max_insert_tries);

					for (pos_t i = 1; i < best_len; ++i)
					{
						if (pos + i + min_length <= (pos_t)raw.size())
						{
							auto h2 = hasher(&raw[pos + i], min_length);
							insert_ht(h2, pos + i, max_offset, max_insert_tries);
						}
					}
					pos += best_len;
				}
				else
				{
					encoder.encode_literal(raw[pos]);
					// Insert the current sequence to the hash table
					if (possible_insert_pos != ~uint32_t(0))
						h = possible_insert_pos;

					insert_ht(h, pos++, max_offset, max_insert_tries);
				}

				if (pos >= next_refresh)
				{
					refresh_ht(raw.data(), pos, max_offset, min_length, max_insert_tries);
					next_refresh = pos + max_offset;
				}
			}

			// Encode any remaining literals
			for (; pos < raw.size(); ++pos)
				encoder.encode_literal(raw[pos]);

			encoder.complete();

			return true;
		}
	};

	// **********************************************************************************
	//
	// **********************************************************************************
	inline bool lzss_rc_compress(std::vector<uint8_t>& raw, std::vector<uint8_t>& packed, uint32_t level, uint8_t max_literal_val = std::numeric_limits<uint8_t>::max())
	{
		refresh::lzss_rc_output ec(packed);

		if (raw.size() > (size_t(1) << 31))
		{
			refresh::lzss_compressor<uint8_t, true> lzss;
			return lzss.compress(raw, ec, level, max_literal_val);
		}
		else
		{
			refresh::lzss_compressor lzss;
			return lzss.compress(raw, ec, level, max_literal_val);
		}
	}
}
