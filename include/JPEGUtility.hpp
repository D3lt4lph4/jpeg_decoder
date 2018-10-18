#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index);
void IDCT(int *new_block, unsigned int component_number);
void YCbCrToBGR(int *new_block);
void FastIDCT(int *new_block, unsigned int component_number);
static void FastIDCT1(int *x, int *y, int ps, int half);

#endif