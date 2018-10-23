#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

#include "JPEGUtils.hpp"

unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index);
void IDCT(int *new_block);
void YCbCrToBGR(JPEGImage *image, std::vector<int> shape);
void FastIDCT(int *new_block);
void FastIDCT1(int *x, int *y, int ps, int half);

#endif