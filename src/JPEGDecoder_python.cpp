#include <pybind11/pybind11.h>

namespace py = pybind11;

void initializeJPEGDecoder(py::module &module) {
  py::module sub_module =
      m.def_submodule("decoder", "The submodule of jpeg_decoder.");

  py::class_<JPEGDecoder> decoder(sub_module, "JPEGDecoder");

  decoder.def("decode_file", &JPEGDecoder::DecodeFile);
}