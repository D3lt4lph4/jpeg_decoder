#include <pybind11/pybind11.h>

namespace py = pybind11;

void initializeJPEGUtils(py::module &);
void initializeJPEGUtility(py::module &);
void initializeJPEGDecoder(py::module &);

PYBIND11_MODULE(jpegdecoder, m) {
  m.doc() = "doc";

  initializeJPEGUtils(m);
  initializeJPEGUtility(m);
  initializeJPEGDecoder(m);
}