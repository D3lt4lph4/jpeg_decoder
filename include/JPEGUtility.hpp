#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

#include "JPEGUtils.hpp"

// Access Functions
unsigned char NextBit(unsigned char *file_content, unsigned int &index,
                      unsigned char &bit_index);

// DCT/IDCT Functions
void IDCT(int *new_block);
void FastIDCT1D(std::vector<int> &x, std::vector<int> &y, const int x_offset,
                const int y_offset, const int ps, const int half,
                const int line_length);
void FastIDCT2D(std::vector<int> &image, const int start_line, const int start_column,
              const int line_length);

// Transformation Functions
void YCbCrToBGR(JPEGImage &image, std::vector<int> shape);

#endif

