#include <emscripten/bind.h>
#include "mdc_reader.h"

using namespace emscripten;
using namespace mdc;


bool query_wrap(query_engine& self, const std::vector<uint32_t>& frame_ids, query_result& res) {
	return self.query(frame_ids, res);
};

//due to some wired reasons I was able to only return this as raw pointer specifically with this type intptr_t
//and on JS side access it via HEAPF32
intptr_t get_coords(frame& self) {
	return (intptr_t)self.coords.data();
}


size_t get_n_coords(const frame& self) {
	return self.coords.size();
}

std::array<std::array<double, 3>, 3> get_box(const frame& self) {
	std::array<std::array<double, 3>, 3> res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res[i][j] = self.box[i][j];
		}
	}
	return res;
}

EMSCRIPTEN_BINDINGS(mdc_wasm) {

	register_vector<std::string>("vec_str");
    register_vector<uint32_t>("vec_u32");
	register_vector<frame>("vec_frame");

	enum_<segment_type>("segment_type")
		.value("unknown", segment_type::unknown)
		.value("molecule", segment_type::molecule)
		.value("other", segment_type::other)
		.value("water", segment_type::water)
		.value("none", segment_type::none);

	value_object<segment_desc_t>("segment_desc_t")
		.field("type", &segment_desc_t::type)
		.field("name", &segment_desc_t::name)
		.field("size", &segment_desc_t::size);


	register_vector<segment_desc_t>("vec_segment_desc_t");

	value_array<std::array<double, 3>>("array_double_3")
        .element(emscripten::index<0>())
        .element(emscripten::index<1>())
		.element(emscripten::index<2>())
        ;

   value_array<std::array<std::array<double, 3>, 3>>("array_array_double_3")
		.element(emscripten::index<0>())
		.element(emscripten::index<1>())
		.element(emscripten::index<2>());

//not needed anymore?
//	class_<atom_coords>("atom_coords")
//		.property("x", &atom_coords::x)
//		.property("y", &atom_coords::y)
//		.property("z", &atom_coords::z)
//		;

	class_<frame>("frame")
		.function("coords", get_coords, allow_raw_pointers())
		.property("n_coords", get_n_coords)
		.property("step", &frame::step)
		.property("time", &frame::time)
		.property("box", get_box)
		.property("prec", &frame::prec)
		;

	class_<query_result>("query_result")
		.constructor<>()
		.property("frames", &query_result::frames)
		;

	class_<reader>("reader")
		.constructor<const std::string&>()
		.function("get_current_error", &reader::get_current_error)
		.function("get_query_engine", &reader::get_query_engine)
		.function("get_no_frames", &reader::get_no_frames)
		.function("get_segments", &reader::get_segments)
		;

	class_<query_engine>("query_engine")
		.function("get_original_atom_ids", &query_engine::get_original_atom_ids)
		.function("query", query_wrap)
		.function("get_current_error", &query_engine::get_current_error)
		;
}
