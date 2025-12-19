#include "engine_decompression.h"
#include <functional>
#include <iomanip>

// *******************************************************************************
void CEngineDecompress::open_stream(const packed_t& packed)
{
	engine_state.clear_rc_dict();
	engine_state.v_rc = packed;
	engine_state.vios->restart_read();
	engine_state.rc->start();
}

// *******************************************************************************
void CEngineDecompress::close_stream()
{
	engine_state.rc->complete();
}

// *******************************************************************************
void CEngineDecompress::decompress(std::shared_ptr<CSegmentModel> model, segment_t& segment, const CSegmentHistory &history)
{
	switch (segment.type)
	{
	case segment_type_t::molecule:
		decompress_molecule(dynamic_pointer_cast<CSegmentModelMolecule>(model), segment, history);
		break;
	case segment_type_t::water:
		decompress_water(dynamic_pointer_cast<CSegmentModelWater>(model), segment, history);
		break;
	case segment_type_t::other:
		decompress_other(dynamic_pointer_cast<CSegmentModelOther>(model), segment, history);
		break;
	case segment_type_t::none:
		// Do nothing here
		break;
	default:
		assert(0);		// !! Never should be here
	}
}

// *******************************************************************************
void CEngineDecompress::decompress_other(std::shared_ptr<CSegmentModelOther> model, segment_t& segment, const CSegmentHistory &history)
{
	segment.coords.clear();

	icoord_t prev_atom;

	std::function<icoord_t(uint32_t i)> predictor;

	if (history.empty())
		predictor = [&](uint32_t i) {return prev_atom; };
	else if (history.size() == 1)
		predictor = [&](uint32_t i) {icoord_t a; history.get(i, a); return history_predictor(a); };
	else if (history.size() == 2)
		predictor = [&](uint32_t i) {icoord_t a, b; history.get(i, a, b); return history_predictor(a, b); };
	else if (history.size() == 3)
		predictor = [&](uint32_t i) {icoord_t a, b, c; history.get(i, a, b, c); return history_predictor(a, b, c); };

	uint32_t ctx_off = history.size() * 10;

	for (uint32_t i = 0; i < segment.n_atoms; ++i)
	{
		auto pred_pos = predictor(i);
		int dx = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dy = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dz = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		icoord_t curr_atom(pred_pos.x + dx, pred_pos.y + dy, pred_pos.z + dz);

		segment.coords.emplace_back(curr_atom);

		prev_atom = curr_atom;
	}
}

// *******************************************************************************
void CEngineDecompress::decompress_water(std::shared_ptr<CSegmentModelWater> model, segment_t& segment, const CSegmentHistory &history)
{
	segment.coords.clear();

	int off_O = model->atom_order;
	int off_H1 = (off_O + 1) % 3;
	int off_H2 = (off_O + 2) % 3;

	auto& coords = segment.coords;
	int64_t d2_OH = model->dist2_OH;
	int64_t d2_HH = model->dist2_HH;
	int64_t d_OH = isqrt(uint64_t(d2_OH));
	int64_t d_HH = isqrt(uint64_t(d2_HH));

	icoord_t atom_O_prev, atom_O;

	std::function<icoord_t(uint32_t i)> predictor_O, predictor_H;

	uint32_t ctx_off = history.size() * 10;

	if (history.empty())
	{
		predictor_O = [&](uint32_t i) {return atom_O_prev; };
		predictor_H = [&](uint32_t i) {return atom_O; };
	}
	else if (history.size() == 1)
	{
		predictor_O = [&](uint32_t i) {icoord_t a; history.get(i, a); return a; };
		predictor_H = [&](uint32_t i) {return atom_O; };
	}
	else if (history.size() == 2)
	{
		predictor_O = [&](uint32_t i) {icoord_t a, b; history.get(i, a, b); return history_predictor(a, b); };
		predictor_H = [&](uint32_t i) {return atom_O; };
	}
	else if (history.size() == 3)
	{
		predictor_O = [&](uint32_t i) {icoord_t a, b, c; history.get(i, a, b, c); return history_predictor(a, b, c); };
		predictor_H = [&](uint32_t i) {return atom_O; };
	}

	for (uint32_t i = 0; i < segment.n_atoms; i += 3)
	{
		auto O_pred_pos = predictor_O(i + off_O);

		// *** O
		int dx_O = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dy_O = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dz_O = rc_decode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		atom_O = icoord_t(O_pred_pos.x + dx_O, O_pred_pos.y + dy_O, O_pred_pos.z + dz_O);
		atom_O_prev = atom_O;

		// *** H1
		auto H1_pred_pos = predictor_H(i + off_H1);
		int dx_H1 = rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		int a = std::clamp((int)(5 * ipow2(H1_pred_pos.x + dx_H1 - atom_O.x) / d2_OH), 0, 5);

//		int dy_H1 = rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dy_H1 = rc_decode_int(ctx_off + 20 + a, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dz_H1 = 0;

		icoord_t atom_H1(H1_pred_pos.x + dx_H1, H1_pred_pos.y + dy_H1, 0);						// z-coord will be set later

		if (compression_features.water_z1 == compression_features_t::water_z1_t::predicted)
		{
			auto pz_H1 = pred_z1(d2_OH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y);
			int dz_H1a = pz_H1.first;
			int dz_H1b = pz_H1.second;

			if (dz_H1a != dz_H1b)
			{
				bool first_predicted = (bool)engine_state.rc->get_cumulative_freq_1bit();
				engine_state.rc->update_frequency(1, (int)first_predicted, 2);

				if (first_predicted)
					dz_H1 = dz_H1a;
				else
					dz_H1 = dz_H1b;
			}
			else
				dz_H1 = dz_H1a;
		}
		else
			dz_H1 = H1_pred_pos.z;

		atom_H1.z = dz_H1 + rc_decode_int(ctx_off + 2, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		// *** H2
		auto H2_pred_pos = predictor_H(i + off_H2);
		if (compression_features.water_x2 == compression_features_t::water_x2_t::predicted)
			H2_pred_pos.x = pred_x2(d_OH, d_HH, d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z);

		int dx_H2 = rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		icoord_t atom_H2(H2_pred_pos.x + dx_H2, 0, 0);									// y- and z-coordinates will be set later

		if (compression_features.water_y2 == compression_features_t::water_y2_t::predicted)
		{
			auto py_H2 = pred_y2_new(d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z, atom_H2.x);
			int dy_H2 = 0;

			if (py_H2.first != py_H2.second)
			{
				bool first_predicted = (bool)engine_state.rc->get_cumulative_freq_1bit();
				engine_state.rc->update_frequency(1, (int)first_predicted, 2);

				if (first_predicted)
					dy_H2 = py_H2.first;
				else
					dy_H2 = py_H2.second;
			}
			else
				dy_H2 = py_H2.first;				

			atom_H2.y = dy_H2 + rc_decode_int(ctx_off + 3, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		}
		else
			atom_H2.y = H2_pred_pos.y + rc_decode_int(ctx_off + 3, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		if (compression_features.water_z2 == compression_features_t::water_z2_t::predicted)
		{
			auto dz_H2 = pred_z2(d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z, atom_H2.x, atom_H2.y);
			atom_H2.z = dz_H2 + rc_decode_int(ctx_off + 2, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		}
		else
			atom_H2.z = H2_pred_pos.z + rc_decode_int(ctx_off + 2, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		uint32_t old_size = coords.size();
		coords.resize(coords.size() + 3);
		coords[old_size + off_O] = atom_O;
		coords[old_size + off_H1] = atom_H1;
		coords[old_size + off_H2] = atom_H2;
	}
}

// *******************************************************************************
void CEngineDecompress::decompress_molecule(std::shared_ptr<CSegmentModelMolecule> model, segment_t& segment, const CSegmentHistory &history)
{
	segment.coords.clear();

	icoord_t prev_atom;

	std::function<icoord_t(uint32_t i)> predictor;

	uint32_t ctx_off = history.size() * 10;

	if (history.empty())
		predictor = [&](uint32_t i) {return prev_atom; };
	else if (history.size() == 1)
		predictor = [&](uint32_t i) {icoord_t a; history.get(i, a); return history_predictor(a); };
	else if (history.size() == 2)
		predictor = [&](uint32_t i) {icoord_t a, b; history.get(i, a, b); return history_predictor(a, b); };
	else if (history.size() == 3)
		predictor = [&](uint32_t i) {icoord_t a, b, c; history.get(i, a, b, c); return history_predictor(a, b, c); };

	for (uint32_t i = 0; i < std::min<uint32_t>(3, segment.n_atoms); ++i)
	{
		auto pred_pos = predictor(i);

		int32_t dx = (int32_t)rc_decode_int(ctx_off + 0, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int32_t dy = (int32_t)rc_decode_int(ctx_off + 0, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int32_t dz = (int32_t)rc_decode_int(ctx_off + 0, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		icoord_t curr_atom(pred_pos.x + dx, pred_pos.y + dy, pred_pos.z + dz);

		segment.coords.emplace_back(curr_atom);

		prev_atom = curr_atom;
	}

	auto& ref_atoms = model->ref_atoms;

//	bool showing = false;

	for (uint32_t i = 3; i < segment.n_atoms; ++i)
	{
		uint32_t ra1 = i - ref_atoms[i].atom_ids[0];
		uint32_t ra2 = i - ref_atoms[i].atom_ids[1];
		uint32_t ra3 = i - ref_atoms[i].atom_ids[2];

		int dx = (int)rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dy = (int)rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);
		int dz = (int)rc_decode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix);

		if (compression_features.molecule == compression_features_t::molecule_t::tetrahedron ||
			(compression_features.molecule == compression_features_t::molecule_t::translate_and_rotate && history.empty()))
		{
			Vec3Dfp c1(segment.coords[ra1].x, segment.coords[ra1].y, segment.coords[ra1].z);
			Vec3Dfp c2(segment.coords[ra2].x, segment.coords[ra2].y, segment.coords[ra2].z);
			Vec3Dfp c3(segment.coords[ra3].x, segment.coords[ra3].y, segment.coords[ra3].z);

			std::pair<Vec3Dfp, Vec3Dfp> predicted;

			if (history.empty())
				predicted = get_q_both(c1, c2, c3, ref_atoms[i].distances[0] * ref_atoms[i].distances[0], ref_atoms[i].distances[1] * ref_atoms[i].distances[1], ref_atoms[i].distances[2] * ref_atoms[i].distances[2], ref_atoms[i].is_left_handed);
			else
				predicted = get_q_both(c1, c2, c3, history.dist2(i, ra1), history.dist2(i, ra2), history.dist2(i, ra3), ref_atoms[i].is_left_handed);

			Vec3Dfp best_prediction;

			int px1 = int_round_or_0(predicted.first.x);
			int py1 = int_round_or_0(predicted.first.y);
			int pz1 = int_round_or_0(predicted.first.z);
			int px2 = int_round_or_0(predicted.second.x);
			int py2 = int_round_or_0(predicted.second.y);
			int pz2 = int_round_or_0(predicted.second.z);

			// !!! Can be ommited if both predictions are equal (after rounding to int)
			if (px1 != px2 || py1 != py2 || pz1 != pz2)
			{
				bool first_predicted;
				if (compression_features.molecule_orientation == compression_features_t::molecule_orientation_t::predicted)
				{
					auto p_binary = rc_find_context(engine_state.dict_binary, 0, *engine_state.tpl_binary);
					first_predicted = p_binary->decode();
				}
				else
				{
					first_predicted = (bool)engine_state.rc->get_cumulative_freq_1bit();
					engine_state.rc->update_frequency(1, (int)first_predicted, 2);
				}

				if (first_predicted)
				{
					dx += px1;
					dy += py1;
					dz += pz1;
				}
				else
				{
					dx += px2;
					dy += py2;
					dz += pz2;
				}
			}
			else
			{
				dx += px1;
				dy += py1;
				dz += pz1;
			}
		}
		else if (compression_features.molecule == compression_features_t::molecule_t::translate_and_rotate)
		{
			Vec3Dfp c1(segment.coords[ra1].x, segment.coords[ra1].y, segment.coords[ra1].z);
			Vec3Dfp c2(segment.coords[ra2].x, segment.coords[ra2].y, segment.coords[ra2].z);
			Vec3Dfp c3(segment.coords[ra3].x, segment.coords[ra3].y, segment.coords[ra3].z);

			icoord_t c1_p, c2_p, c3_p, q_p;
			history.get(ra1, c1_p);
			history.get(ra2, c2_p);
			history.get(ra3, c3_p);
			history.get(i, q_p);
			Vec3Dfp c1_prev(c1_p.x, c1_p.y, c1_p.z), c2_prev(c2_p.x, c2_p.y, c2_p.z), c3_prev(c3_p.x, c3_p.y, c3_p.z), q_prev(q_p.x, q_p.y, q_p.z);

			Vec3Dfp predicted = get_q_new2(c2_prev, c3_prev, c1_prev, c2, c3, c1, q_prev);

			if (predicted.is_wrong())
			{
				predicted.x = q_p.x + segment.coords[ra1].x - c1_p.x;
				predicted.y = q_p.y + segment.coords[ra1].y - c1_p.y;
				predicted.z = q_p.z + segment.coords[ra1].z - c1_p.z;
			}

			dx += int_round_or_0(predicted.x);
			dy += int_round_or_0(predicted.y);
			dz += int_round_or_0(predicted.z);
		}
		else if (compression_features.molecule == compression_features_t::molecule_t::moving_3vectors && !history.empty())
		{
			icoord_t p_i, c_i;
			icoord_t p_ra1, p_ra2, p_ra3;
			icoord_t c_ra1, c_ra2, c_ra3;

			history.get(i, p_i);
			history.get(ra1, p_ra1);
			history.get(ra2, p_ra2);
			history.get(ra3, p_ra3);
			c_ra1 = segment.coords[ra1];
			c_ra2 = segment.coords[ra2];
			c_ra3 = segment.coords[ra3];

			c_i.x = p_i.x + (12 * (c_ra1.x - p_ra1.x) + 3 * (c_ra2.x - p_ra2.x) + (c_ra3.x - p_ra3.x)) / 16;
			c_i.y = p_i.y + (12 * (c_ra1.y - p_ra1.y) + 3 * (c_ra2.y - p_ra2.y) + (c_ra3.y - p_ra3.y)) / 16;
			c_i.z = p_i.z + (12 * (c_ra1.z - p_ra1.z) + 3 * (c_ra2.z - p_ra2.z) + (c_ra3.z - p_ra3.z)) / 16;

			dx += c_i.x;
			dy += c_i.y;
			dz += c_i.z;
		}
		else if (compression_features.molecule == compression_features_t::molecule_t::moving_1vector && !history.empty())
		{
			icoord_t p_i, c_i;
			icoord_t p_ra1;
			icoord_t c_ra1;

			history.get(i, p_i);
			history.get(ra1, p_ra1);
			c_ra1 = segment.coords[ra1];

			c_i.x = p_i.x + c_ra1.x - p_ra1.x;
			c_i.y = p_i.y + c_ra1.y - p_ra1.y;
			c_i.z = p_i.z + c_ra1.z - p_ra1.z;

			dx += c_i.x;
			dy += c_i.y;
			dz += c_i.z;
		}
		else
		{
			auto predicted = predictor(i);

			dx += predicted.x;
			dy += predicted.y;
			dz += predicted.z;
		}

		segment.coords.emplace_back(icoord_t(dx, dy, dz));
	}
}

// *******************************************************************************
void CEngineDecompress::decompress_raw(segment_t& segment)
{
	auto& packed = segment.packed;
	auto& coords = segment.coords;

	uint32_t n_coords = (uint32_t) packed.size() / 12;

	coords.resize(n_coords);

	auto p = packed.begin();

	for (uint32_t i = 0; i < n_coords; ++i)
	{
		for (int& a : { std::ref(coords[i].x), std::ref(coords[i].y), std::ref(coords[i].z) })
		{
			uint32_t una = 0;

			for (int j = 0; j < 4; ++j)
				una += ((uint32_t)(*p++)) << (8 * j);

			a = (int)una;
		}
	}

	segment.n_atoms = n_coords;

	packed.clear();
	packed.shrink_to_fit();
}

// *****************************************************************
int64_t CEngineDecompress::rc_decode_uint(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix)
{
//	bool is_neg = false;
	int val = 0;
	int nb;

	if (compression_features.rc_uint32 == compression_features_t::rc_uint32_t::prefix_fast)
	{
		nb = (int)engine_state.rc->get_cumulative_freq(engine_state.prefix_max_bits + 2);
		engine_state.rc->update_frequency(1, nb, engine_state.prefix_max_bits + 2);
	}
	else
	{
		auto p_prefix = rc_find_context(dict_prefix, ctx, *engine_state.tpl_prefix);
		nb = (int) p_prefix->decode();
	}

	if (nb <= engine_state.prefix_max_bits)
	{
		if (nb <= 1)
			return nb;

		if (compression_features.rc_uint32 == compression_features_t::rc_uint32_t::suffix_exact && nb < engine_state.prefix_short_max_bits)
		{
			auto p_suffix = rc_find_context(dict_suffix, nb, &engine_state.tpl_suffix[nb - 1]);
			val = p_suffix->decode();
		}
		else
		{
			val = (int)engine_state.rc->get_cumulative_freq_pow2(nb - 1);
			engine_state.rc->update_frequency(1, val, 1ull << (nb - 1));
		}

		val += (1 << (nb - 1));
	}
	else
	{
		auto big_context = rc_find_context(dict_suffix, ctx * 32 + 31, &engine_state.tpl_suffix.back());

		uint32_t b1 = rc_find_context(dict_suffix, ctx * 32 + 29, &engine_state.tpl_suffix.back())->decode();
		uint32_t b2 = rc_find_context(dict_suffix, ctx * 32 + 30, &engine_state.tpl_suffix.back())->decode();
		uint32_t b3 = big_context->decode();
		uint32_t b4 = 0;

		if (b3 == 255)
		{
			b3 = big_context->decode();
			b4 = big_context->decode();
		}

		val = (b4 << 24) + (b3 << 16) + (b2 << 8) + b1;

		val += 1 << (engine_state.prefix_max_bits - 1);
	}

	return val;
}

// *****************************************************************
int64_t CEngineDecompress::rc_decode_int(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix)
{
	int64_t val = rc_decode_uint(ctx, dict_prefix, dict_suffix);

	if (val & 1)
		val = -(val + 1) / 2;
	else
		val /= 2;

	return val;
}
