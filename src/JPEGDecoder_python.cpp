#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "JPEGDecoder.hpp"
#include "JPEGUtils.hpp"

namespace py = pybind11;

void initializeJPEGDecoder(py::module &module) {
  py::module submodule = module.def_submodule("decoder");

  py::class_<JPEGDecoder>(submodule, "JPEGDecoder")
      .def(py::init<>())
      .def(py::init<const unsigned char>())
      .def("decode_file", &JPEGDecoder::DecodeFile,
           py::return_value_policy::reference)
      .def("get_image_size_x", &JPEGDecoder::getImageSizeX)
      .def("get_image_size_y", &JPEGDecoder::getImageSizeY)
      .def("get_block_per_line", &JPEGDecoder::getBlockPerLine)
      .def("get_block_per_column", &JPEGDecoder::getBlockPerColumn)
      .def("get_channels", &JPEGDecoder::getChannels);
}
