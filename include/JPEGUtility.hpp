#ifndef __JPEG_UTILITY__
#define __JPEG_UTILITY__

unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index);
void IDCT(int *new_block);
void YCbCrToBGR(int *first_component, int *second_component, int *third_component);
void FastIDCT(int *new_block);
void FastIDCT1(int *x, int *y, int ps, int half);

#endif