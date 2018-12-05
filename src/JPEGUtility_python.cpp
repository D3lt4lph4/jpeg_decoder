/**
 * \file JPEGUtility.cpp
 * \author Deguerre Benjamin
 * \brief Contains all the utility functions for the JPEGDecoder.
 */

#include "JPEGUtility.hpp"



#include <math.h>
#include <iostream>
#include <stdexcept>

#include "JPEGType.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;


void initializeJPEGUtility(py::module &module) {
  py::module sub_module = module.def_submodule("utility", "The submodule of jpeg_decoder.");

  sub_module.def("next_bit", &NextBit);


  sub_module.def("fast_idct_1d", &FastIDCT1D);

  sub_module.def("fast_idct_2d", &FastIDCT2D);

  sub_module.def("fast_idct_row", &IDCTRow);

  sub_module.def("fast_idct_col", &IDCTCol);

  sub_module.def("fast_idct_2d_second", &FastIDCT2D_Second);

  sub_module.def("de_level_shift", &DeLevelShift);

  sub_module.def("YCbCr_to_RGB", &YCbCrToBGR);

}
