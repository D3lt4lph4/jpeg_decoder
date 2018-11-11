#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

#include "JPEGUtils.hpp"

#define W1 2841 /* 2048*sqrt(2)*cos(1*pi/16) */
#define W2 2676 /* 2048*sqrt(2)*cos(2*pi/16) */
#define W3 2408 /* 2048*sqrt(2)*cos(3*pi/16) */
#define W5 1609 /* 2048*sqrt(2)*cos(5*pi/16) */
#define W6 1108 /* 2048*sqrt(2)*cos(6*pi/16) */
#define W7 565 /* 2048*sqrt(2)*cos(7*pi/16) */

// Access Functions
unsigned char NextBit(unsigned char *file_content, unsigned int &index,
                      unsigned char &bit_index);

// DCT/IDCT Functions
void IDCT(int *new_block);
void FastIDCT1D(std::vector<int> &x, std::vector<int> &y, const int x_offset,
                const int y_offset, const int ps, const int half,
                const int line_length);
void FastIDCT2D(std::vector<int> &image, const int start_line,
                const int start_column, const int line_length);

// Other implementation of the fast IDCT 2D
void IDCT_Row(std::vector<int> &x, const int x_offset, const int y_offset,
              const int line_length);
void IDCT_Col(std::vector<int> &x, const int x_offset, const int y_offset,
              const int line_length);
void FastIDCT2D_Second(std::vector<int> &image, const int start_line,
                       const int start_column, const int line_length);

// Transformation Functions
void YCbCrToBGR(JPEGImage &image, std::vector<int> shape);
void DeLevelShift(JPEGImage &image);
#endif
