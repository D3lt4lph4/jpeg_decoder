#include "JPEGUtility.hpp"

#include <math.h>
#include <iostream>
#include <stdexcept>

#include "JPEGType.hpp"

#define xadd3(xa, xb, xc, xd, h)                                               \
  p = xa + xb, n = xa - xb, xa = p + xc + h, xb = n + xd + h, xc = p - xc + h, \
  xd = n - xd + h  // triple-butterfly-add (and possible rounding)
#define xmul(xa, xb, k1, k2, sh)                               \
  n = k1 * (xa + xb), p = xa, xa = (n + (k2 - k1) * xb) >> sh, \
  xb = (n - (k2 + k1) * p) >> sh  // butterfly-mul equ.(2)


/**
 * \fn unsigned char NextBit()
 * \brief Returns the next bit in the stream.
 *
 * This function extract the bit at the specified position in the stream. When
 * the bit_index value reach 0, the bit_index is set to 8 and the index is
 * incremented by 1. If a value outside of [1 ; 8] is provided for the
 * bit_index, the function will throw an out_of_range error.
 *
 * If a 0xFF byte if found in the stream, if the following byte is 0x00, the bit
 * is extracted, else an error is thrown. If 0xFF00 is encountered, when
 * accessing the lowest bit in the 0xFF, the index will be incremented by 2
 * instead of 1.
 *
 * \param[in] file_content Pointer to the unsigned char to retrieve the bit
 * from. \param[in, out] index Pointer containing the offset for the pointer to
 * the byte to extract the bit from. \param[in, out] bit_index The index of the
 * bit to extract. A value of 8 extract the highest bit in the byte, a value of
 * 1 extract the lowest.
 *
 * \return Either 0 or 1, the value of the bit.
 */
unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index) {
  unsigned char current_byte = file_content[*index], bit;
  if (*bit_index < 1 || *bit_index > 8) {
    throw std::out_of_range(
        "The index provided for the bit is out of range with a value of " +
        *bit_index);
  }
  if (*bit_index == 8) {
    if (current_byte == 0xFF) {
      if (!(file_content[*index + 1] == 0x00)) {
        if (file_content[*index + 1] == DEFINE_NUMBER_OF_LINE) {
          throw std::runtime_error("DNL byte found");
        } else {
          throw std::runtime_error(
              "Error while parsing the file, DNL byte expected but something "
              "else found");
        }
      }
    }
  }

  bit = current_byte << (8 - *bit_index);
  bit = bit >> 7;
  *bit_index -= 1;

  if (*bit_index == 0) {
    if (current_byte == 0xFF) {
      *index += 2;
    } else {
      *index += 1;
    }
    *bit_index = 8;
  }
  return bit;
}

/**
 * \fn void IDCT(cv::Mat *new_block, unsigned int component_number)
 *
 * \brief Calculate the inverse cosinus transform over the block provided. The
 * transform is calculated the channel provided in component_number.
 *
 * \param[in, out] new_block A pointer to the OpenCV data to modify.
 *
 * \param[in] component_number The component on which to perform the IDCT.
 */
void IDCT(int *new_block) {
  float result, cu, cv;
  int temp_operation[64];

  
  for(size_t c = 0; c < 8; c++)
  {
    for(size_t r = 0; r < 8; r++)
    {
      temp_operation[c*8+r] = new_block[c*8+r];
    }
  }

  for (size_t x = 0; x < 8; x++) {
    for (size_t y = 0; y < 8; y++) {
      result = 0;
      for (size_t u = 0; u < 8; u++) {
        for (size_t v = 0; v < 8; v++) {
          if (u == 0) {
            cu = 1 / sqrt(2);
          } else {
            cu = 1;
          }
          if (v == 0) {
            cv = 1 / sqrt(2);
          } else {
            cv = 1;
          }
          result += cu * cv *
                    temp_operation[v+u*8] *
                    cos((2 * x + 1) * u * M_PI / 16.0) *
                    cos((2 * y + 1) * v * M_PI / 16.0);
        }
      }
      result = result / 4 + 128;
      if (result > 255) {
        new_block[y+x*8] = 255;
      } else if (result < 0) {
        new_block[y+x*8] = 0;
      } else {
        new_block[y+x*8] = (int)result;
      }
    }
  }
}

/**
 * \fn void FastIDCT1(int *x, int *y, int ps, int half)
 * \brief Compute the one dimension IDCT.
 * 
 * \param[in,out] x, no se
 * \param[in,out] y, no se
 * \param[in] ps, no se
 * \param[in] half, no se
 */
void FastIDCT1(int *x, int *y, int ps, int half)
{
  int p, n;
  x[0] <<= 9, x[1] <<= 7, x[3] *= 181, x[4] <<= 9, x[5] *= 181, x[7] <<= 7;
  xmul(x[6], x[2], 277, 669, 0);
  xadd3(x[0], x[4], x[6], x[2], half);
  xadd3(x[1], x[7], x[3], x[5], 0);
  xmul(x[5], x[3], 251, 50, 6);
  xmul(x[1], x[7], 213, 142, 6);
  y[0 * 8] = (x[0] + x[1]) >> ps;
  y[1 * 8] = (x[4] + x[5]) >> ps;
  y[2 * 8] = (x[2] + x[3]) >> ps;
  y[3 * 8] = (x[6] + x[7]) >> ps;
  y[4 * 8] = (x[6] - x[7]) >> ps;
  y[5 * 8] = (x[2] - x[3]) >> ps;
  y[6 * 8] = (x[4] - x[5]) >> ps;
  y[7 * 8] = (x[0] - x[1]) >> ps;
}

/**
 * \fn void FastIDCT(int *new_block, unsigned int component_number)
 * \brief Compute the a faster IDCT in dimension two. The algorithm is based
 * upon the DCT implementation of Loeffler.
 *
 * \param[in] int* new_block, An array of size 8x8 representing the coefficients
 * of the DCT of a block.
 *
 *
 */
void FastIDCT(int *new_block)  // 2D 8x8 IDCT
{
  int i, b[64], b2[64];

  for (i = 0; i < 8; i++)
    FastIDCT1(new_block + i * 8, b2 + i, 9, 1 << 8);  // row
  for (i = 0; i < 8; i++)
    FastIDCT1(b2 + i * 8, new_block + i, 12, 1 << 11);  // col

  for (size_t c = 0; c < 8; c++) {
    for (size_t r = 0; r < 8; r++) {
      new_block[c * 8 + r] = new_block[c * 8 + r] + 128;
      if (new_block[c * 8 + r] > 255) {
        new_block[c * 8 + r] = 255;
      } else if (new_block[c * 8 + r] < 0) {
        new_block[c * 8 + r] = 0;
      } else {
        new_block[c * 8 + r] = (int)new_block[c * 8 + r];
      }
    }
  }
}

/**
 * \fn void YCbCrToBGR(int *new_block)
 * \brief This function transform a "block" of size 8*8*3 from YCbCr space to RGB space.
 * 
 * \param[in,out] int *new_block, a pointer to the block of data which will be subject the transformation
 * 
 */
void YCbCrToBGR(int *first_component, int *second_component, int *third_component) {
  int R, G, B;
  for (size_t row = 0; row < 8; row++) {
    for (size_t col = 0; col < 8; col++) {
      R = first_component[row*8+col] + 1.402 * (third_component[row*8+col+128] - 128);
      G = first_component[row*8+col] - 0.34414 * (second_component[row*8+col+64] - 128) - 0.71414 * (third_component[row*8+col+128] - 128);
      B = first_component[row*8+col] + 1.772 * (second_component[row*8+col+64] - 128);
      third_component[row*8+col+128] = R;
      second_component[row*8+col+64] = G;
      first_component[row*8+col] = B;
    }
  }
}