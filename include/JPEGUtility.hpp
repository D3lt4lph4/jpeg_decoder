#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

#include <opencv2/core/core.hpp>

unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index);
void IDCT(cv::Mat *new_block, unsigned int component_number);
void YCbCrToBGR(cv::Mat *new_block);
#endif