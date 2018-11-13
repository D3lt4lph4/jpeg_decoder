#include "JPEGDecoder.hpp"
#include "JPEGUtils.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

void initializeJPEGDecoder(py::module &module) {
  py::module sub_module =
      module.def_submodule("decoder", "The submodule of jpeg_decoder.");

  py::class_<JPEGDecoder>(sub_module, "JPEGDecoder")
      .def(py::init<>())
      .def(py::init<const unsigned char>())
      .def("decode_file", &JPEGDecoder::DecodeFile, py::return_value_policy::reference)
      .def("get_image_size_x", &JPEGDecoder::getImageSizeX)
      .def("get_image_size_y", &JPEGDecoder::getImageSizeY);
}
