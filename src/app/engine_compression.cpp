// compress - i trajectory.xtc - o trajectory.mdc - d trajectory.desc - a 10

#include "engine_compression.h"
#include <functional>
#include <algorithm>

// *******************************************************************************
void CEngineCompress::open_stream()
{
	engine_state.rc->start();
}

// *******************************************************************************
void CEngineCompress::close_stream(packed_t& packed)
{
	engine_state.rc->complete();
	packed = engine_state.v_rc;
	engine_state.clear_rc_dict();
}

// *******************************************************************************
void CEngineCompress::compress(const std::shared_ptr<CSegmentModel> model, const segment_t& segment, const CSegmentHistory &history)
{
	switch (segment.type)
	{
	case segment_type_t::molecule:
		compress_molecule(dynamic_pointer_cast<CSegmentModelMolecule>(model), segment, history);
		break;
	case segment_type_t::water:
		compress_water(dynamic_pointer_cast<CSegmentModelWater>(model), segment, history);
		break;
	case segment_type_t::other:
		compress_other(dynamic_pointer_cast<CSegmentModelOther>(model), segment, history);
		break;
	case segment_type_t::none:
		// Do nothing here
		break;
	default:
		assert(0);		// !! Never should be here
	}
}

// *******************************************************************************
void CEngineCompress::compress_other(std::shared_ptr<CSegmentModelOther> model, const segment_t& segment, const CSegmentHistory &history)
{
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

	for(uint32_t i = 0; i < segment.coords.size(); ++i)
	{
		auto& atom = segment.coords[i];
		auto pred_pos = predictor(i);
		int dx = atom.x - pred_pos.x;
		int dy = atom.y - pred_pos.y;
		int dz = atom.z - pred_pos.z;

		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz);

		prev_atom = atom;
	}
}

// *******************************************************************************
void CEngineCompress::compress_water(std::shared_ptr<CSegmentModelWater> model, const segment_t& segment, const CSegmentHistory &history)
{
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
		atom_O = coords[i + off_O];
		auto& atom_H1 = coords[i + off_H1];
		auto& atom_H2 = coords[i + off_H2];

		auto O_pred_pos = predictor_O(i + off_O);	

		// *** O
		int dx_O = atom_O.x - O_pred_pos.x;
		int dy_O = atom_O.y - O_pred_pos.y;
		int dz_O = atom_O.z - O_pred_pos.z;

		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx_O);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy_O);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz_O);
		atom_O_prev = atom_O;

		// *** H1
		auto H1_pred_pos = predictor_H(i + off_H1);
		int dx_H1 = atom_H1.x - H1_pred_pos.x;
		int dy_H1 = atom_H1.y - H1_pred_pos.y;
		int dz_H1 = 0;

		int a = std::clamp((int) (5 * ipow2(atom_H1.x - atom_O.x) / d2_OH), 0, 5);

		rc_encode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx_H1);
		rc_encode_int(ctx_off + 20 + a, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy_H1);

		if (compression_features.water_z1 == compression_features_t::water_z1_t::predicted)
		{
			auto pz_H1 = pred_z1(d2_OH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y);
			int dz_H1a = atom_H1.z - pz_H1.first;
			int dz_H1b = atom_H1.z - pz_H1.second;

			if (dz_H1a != dz_H1b)
			{
				if (iabs(dz_H1a) < iabs(dz_H1b))
				{
					engine_state.rc->encode_1bit(1);
					dz_H1 = dz_H1a;
				}
				else
				{
					engine_state.rc->encode_1bit(0);
					dz_H1 = dz_H1b;
				}
			}
			else
				dz_H1 = dz_H1a;
		}
		else
			dz_H1 = atom_H1.z - H1_pred_pos.z;

		rc_encode_int(ctx_off + 2, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz_H1);

		// *** H2
		auto H2_pred_pos = predictor_H(i + off_H2);

		if (compression_features.water_x2 == compression_features_t::water_x2_t::predicted)
			H2_pred_pos.x = pred_x2(d_OH, d_HH, d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z);

		int dx_H2 = atom_H2.x - H2_pred_pos.x;

		rc_encode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx_H2);

		std::pair<int, int> py_H2{};
		if (compression_features.water_y2 == compression_features_t::water_y2_t::predicted)
			py_H2 = pred_y2_new(d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z, atom_H2.x);
		else
			py_H2 = std::make_pair(H2_pred_pos.y, H2_pred_pos.y);

		int dy_H2a = atom_H2.y - py_H2.first;
		int dy_H2b = atom_H2.y - py_H2.second;
		int dy_H2 = 0;

		if (py_H2.first != py_H2.second)
		{
			if (iabs(dy_H2a) < iabs(dy_H2b))
			{
				engine_state.rc->encode_1bit(1);
				dy_H2 = dy_H2a;
			}
			else
			{
				engine_state.rc->encode_1bit(0);
				dy_H2 = dy_H2b;
			}
		}
		else
		{
			dy_H2 = dy_H2a;
		}

		rc_encode_int(ctx_off + 3, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy_H2);

		int pz_H2 = 0;
		if (compression_features.water_z2 == compression_features_t::water_z2_t::predicted)
			pz_H2 =	pred_z2(d2_OH, d2_HH, atom_O.x, atom_O.y, atom_O.z, atom_H1.x, atom_H1.y, atom_H1.z, atom_H2.x, atom_H2.y);
		else
			pz_H2 = H2_pred_pos.z;
		int dz_H2 = atom_H2.z - pz_H2;

		rc_encode_int(ctx_off + 2, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz_H2);
	}
}

// *******************************************************************************
void CEngineCompress::compress_molecule(std::shared_ptr<CSegmentModelMolecule> model, const segment_t& segment, const CSegmentHistory &history)
{
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
		const icoord_t &atom = segment.coords[i];
		auto pred_pos = predictor(i);

		int dx = atom.x - pred_pos.x;
		int dy = atom.y - pred_pos.y;
		int dz = atom.z - pred_pos.z;

		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy);
		rc_encode_int(ctx_off, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz);

		prev_atom = atom;
	}

	auto& ref_atoms = model->ref_atoms;

	for (uint32_t i = 3; i < segment.coords.size(); ++i)
	{
		uint32_t ra1 = i - ref_atoms[i].atom_ids[0];
		uint32_t ra2 = i - ref_atoms[i].atom_ids[1];
		uint32_t ra3 = i - ref_atoms[i].atom_ids[2];
		int bit_to_encode = -1;

		int dx, dy, dz;

		if (compression_features.molecule == compression_features_t::molecule_t::tetrahedron || 
			(compression_features.molecule == compression_features_t::molecule_t::translate_and_rotate && history.empty()))
		{
			Vec3Dfp c1(segment.coords[ra1].x, segment.coords[ra1].y, segment.coords[ra1].z);
			Vec3Dfp c2(segment.coords[ra2].x, segment.coords[ra2].y, segment.coords[ra2].z);
			Vec3Dfp c3(segment.coords[ra3].x, segment.coords[ra3].y, segment.coords[ra3].z);
			Vec3Dfp q(segment.coords[i].x, segment.coords[i].y, segment.coords[i].z);

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

			dx = segment.coords[i].x;
			dy = segment.coords[i].y;
			dz = segment.coords[i].z;
			bool first;

			if (vector_length2(q - predicted.first) < vector_length2(q - predicted.second))
			{
				dx -= px1;
				dy -= py1;
				dz -= pz1;
				first = true;
			}
			else
			{
				dx -= px2;
				dy -= py2;
				dz -= pz2;
				first = false;
			}

			// !!! Can be ommited if both predictions are equal (after rounding to int)
			if (px1 != px2 || py1 != py2 || pz1 != pz2)
				bit_to_encode = (int) first;
		}
		else if (compression_features.molecule == compression_features_t::molecule_t::translate_and_rotate)
		{
			Vec3Dfp c1(segment.coords[ra1].x, segment.coords[ra1].y, segment.coords[ra1].z);
			Vec3Dfp c2(segment.coords[ra2].x, segment.coords[ra2].y, segment.coords[ra2].z);
			Vec3Dfp c3(segment.coords[ra3].x, segment.coords[ra3].y, segment.coords[ra3].z);
			Vec3Dfp q(segment.coords[i].x, segment.coords[i].y, segment.coords[i].z);

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

//			dx = segment.coords[i].x - iround<int>(predicted.x);
//			dy = segment.coords[i].y - iround<int>(predicted.y);
//			dz = segment.coords[i].z - iround<int>(predicted.z);
			dx = segment.coords[i].x - int_round_or_0(predicted.x);
			dy = segment.coords[i].y - int_round_or_0(predicted.y);
			dz = segment.coords[i].z - int_round_or_0(predicted.z);
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

			dx = segment.coords[i].x - c_i.x;
			dy = segment.coords[i].y - c_i.y;
			dz = segment.coords[i].z - c_i.z;
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

			dx = segment.coords[i].x - c_i.x;
			dy = segment.coords[i].y - c_i.y;
			dz = segment.coords[i].z - c_i.z;
		}
		else
		{
			auto predicted = predictor(i);

			dx = segment.coords[i].x - predicted.x;
			dy = segment.coords[i].y - predicted.y;
			dz = segment.coords[i].z - predicted.z;
		}

		rc_encode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dx);
		rc_encode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dy);
		rc_encode_int(ctx_off + 1, engine_state.dict_cart_prefix, engine_state.dict_cart_suffix, dz);

		if (bit_to_encode != -1)
		{
			if (compression_features.molecule_orientation == compression_features_t::molecule_orientation_t::predicted)
			{
				auto p_binary = rc_find_context(engine_state.dict_binary, 0, *engine_state.tpl_binary);
				p_binary->encode((uint32_t)bit_to_encode);
			}
			else
				engine_state.rc->encode_1bit(bit_to_encode);
		}
	}
}

// *******************************************************************************
void CEngineCompress::rc_encode_int(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix, int64_t val)
{
	val *= 2;
	if (val < 0)
		val = -val - 1;

	rc_encode_uint(ctx, dict_prefix, dict_suffix, (uint64_t)val);
}

// *******************************************************************************
void CEngineCompress::rc_encode_uint(uint32_t ctx, refresh::rc_context_vec_emb<model_prefix_t>& dict_prefix, refresh::rc_context_vec<model_suffix_t>& dict_suffix, uint64_t val)
{
	auto nb = no_bits(val);

	if (nb <= engine_state.prefix_max_bits)
	{
		if (compression_features.rc_uint32 == compression_features_t::rc_uint32_t::prefix_fast)
			engine_state.rc->encode_frequency(1, nb, engine_state.prefix_max_bits + 2);
		else
		{
			auto p_prefix = rc_find_context(dict_prefix, ctx, *engine_state.tpl_prefix);
			p_prefix->encode(nb);
		}

		if (nb <= 1)
			return;

		val -= (((int64_t)1) << (nb - 1));

		if (compression_features.rc_uint32 == compression_features_t::rc_uint32_t::suffix_exact && nb < engine_state.prefix_short_max_bits)
		{
			auto p_suffix = rc_find_context(dict_suffix, nb, &engine_state.tpl_suffix[nb - 1]);
			p_suffix->encode((uint32_t)val);
		}
		else
			engine_state.rc->encode_frequency_pow2(1, val, nb - 1);
	}
	else
	{
		if (compression_features.rc_uint32 == compression_features_t::rc_uint32_t::prefix_fast)
			engine_state.rc->encode_frequency(1, engine_state.prefix_max_bits + 1, engine_state.prefix_max_bits + 2);
		else
		{
			auto p_prefix = rc_find_context(dict_prefix, ctx, *engine_state.tpl_prefix);
			p_prefix->encode(engine_state.prefix_max_bits + 1);
		}

		val -= 1ull << (engine_state.prefix_max_bits - 1);

		int b1 = (int)(val & 0xff);
		int b2 = (int)((val >> 8) & 0xff);
		int b3 = (int)((val >> 16) & 0xff);
		int b4 = (int)(val >> 24);

		rc_find_context(dict_suffix, ctx * 32 + 29, &engine_state.tpl_suffix.back())->encode(b1);
		rc_find_context(dict_suffix, ctx * 32 + 30, &engine_state.tpl_suffix.back())->encode(b2);

		if (b4 == 0 && b3 < 255)
			rc_find_context(dict_suffix, ctx * 32 + 31, &engine_state.tpl_suffix.back())->encode(b3);
		else
		{
			auto big_context = rc_find_context(dict_suffix, ctx * 32 + 31, &engine_state.tpl_suffix.back());

			big_context->encode(255);
			big_context->encode(b3);
			big_context->encode(b4);
		}
	}
}

