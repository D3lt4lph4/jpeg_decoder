#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "JPEGUtils.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace py = pybind11;


void initializeJPEGUtils(py::module &module) {
  py::module sub_module =
      module.def_submodule("utils", "The submodule of jpeg_decoder.");

  py::bind_vector<std::vector<int>>(module, "intVector");

  py::class_<JPEGImage>(sub_module, "JPEGImage")
      .def(py::init<>())
      .def(py::init<std::vector<std::pair<int, int>>>())
      .def("get_component_shape", &JPEGImage::GetComponentShape)
      .def("get_real_shape", &JPEGImage::GetRealShape)
      .def("get_number_of_component", &JPEGImage::GetNumberOfComponent)
      .def("set_real_shape", &JPEGImage::SetRealShape)
      .def("rescale_to_real_size", &JPEGImage::RescaleToRealSize)
      .def("at", &JPEGImage::at, py::return_value_policy::reference)
      .def("get_data", [](JPEGImage &image, int component) {
          auto data = &image.GetData(component);
          return py::array(data->size(), data->data());
      });
}
