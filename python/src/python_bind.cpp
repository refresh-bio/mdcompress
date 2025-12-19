#include "../../src/mdc_lib/mdc_reader.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

NB_MODULE(mdc, m)
{
    nb::enum_<mdc::segment_type>(m, "SegmentType")
        .value("unknown", mdc::segment_type::unknown)
        .value("molecule", mdc::segment_type::molecule)
        .value("other", mdc::segment_type::other)
        .value("water", mdc::segment_type::water)
        .value("none", mdc::segment_type::none)
    ;

    nb::class_<mdc::segment_desc_t>(m, "SegmentDesc")
        .def_ro("type", &mdc::segment_desc_t::type)
        .def_ro("name", &mdc::segment_desc_t::name)
        .def_ro("size", &mdc::segment_desc_t::size)
    ;

    nb::class_<mdc::frame>(m, "Frame")
        .def_ro("step", &mdc::frame::step)
        .def_ro("time", &mdc::frame::time)
        .def_prop_ro("box", [](const mdc::frame& self) {
            return nb::ndarray<float, nb::numpy, nb::shape<3, 3>>((float*)self.box);
        })
        .def_ro("prec", &mdc::frame::prec)
        .def_prop_ro("coords", [](const mdc::frame &self) {
            //if coord type ever change just file compilation because we here cast it to raw float*
            static_assert(3 * sizeof(float) == sizeof(self.coords[0]));

            return nb::ndarray<float, nb::numpy, nb::shape< -1, 3>>(
                (float*)self.coords.data(),
                { self.coords.size(), 3 }
            );
        })
    ;

    nb::class_<mdc::query_result>(m, "QueryResult")
        .def(nb::init<>())
        .def_ro("frames", &mdc::query_result::frames)
    ;

    nb::class_<mdc::query_engine>(m, "QueryEngine")
        .def("get_original_atom_ids", &mdc::query_engine::get_original_atom_ids, nb::rv_policy::reference)
        .def("query", [](mdc::query_engine &self, nb::ndarray<uint32_t, nb::shape<-1>> py_frame_ids, mdc::query_result &result)
             { return self.query(std::span<const uint32_t>(py_frame_ids.data(), py_frame_ids.size()), result); })

        .def("get_current_error", &mdc::query_engine::get_current_error)
    ;

    nb::class_<mdc::reader>(m, "Reader")
        .def(nb::init<const std::string &>())
        .def("get_current_error", &mdc::reader::get_current_error)
        .def("get_segments", &mdc::reader::get_segments, nb::rv_policy::reference)
        .def("get_no_frames", &mdc::reader::get_no_frames)
        .def("get_anchor_ids", &mdc::reader::get_anchor_ids)
        .def("get_query_engine", &mdc::reader::get_query_engine)
    ;

}

