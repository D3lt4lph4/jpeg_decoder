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

/**
 * \def xadd3(xa, xb, xc, xd, h)
 *
 * \brief A macro to realise the set of additions as described <a
 * href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.205.9199&rep=rep1&type=pdf">here</a>
 * (figure 3)
 *
 * The operation is a "triple butterfly addition", it is used twice in the
 * figure 3 in the link above, with two inputs X1 and X2, a simple butterfly
 * addition is X1 = X1 + X2 & X2 = X1 - X2.
 *
 *  This macro changes the input values to:
 *  <ul>
 *    <li>xa = xa + xb + xc</li>
 *    <li>xb = xa - xb + xd</li>
 *    <li>xc = xa + xb - xc</li>
 *    <li>xd = xa - xb - xd</li>
 *  </ul>
 *
 * \param[in, out] xa In figure 3 of the link, would be either X0 or X1.
 * \param[in, out] xb In figure 3 of the link, would be either X4 or X7.
 * \param[in, out] xc In figure 3 of the link, would be either X6 or X3.
 * \param[in, out] xd In figure 3 of the link, would be either X2 or X5.
 */
#define xadd3(xa, xb, xc, xd) \
  p = xa + xb, n = xa - xb, xa = p + xc, xb = n + xd, xc = p - xc, xd = n - xd

/**
 * \def xmul(xa, xb, k1, k2, sh)
 *
 * \brief A macro to realise the set of multiplications as described <a
 * href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.205.9199&rep=rep1&type=pdf">here
 * </a> (figure 3)
 *
 * The operation is a "butterfly multiplication" it is used 3 times in the
 * algorithm used. This macro changes the input values to: <ul> <li>xa = alpha *
 * xa - beta * xb = -(alpha + beta) * xb + alpha * (xa + xb) </li> <li>xb =
 * alpha * xb + beta * xa = (beta - alpha) * xa + alpha * (xa + xb)</li>
 *  </ul>
 *
 * The parameter sh is a parameter to shift the values back to some results to
 * be able to carry out the additions. It could be moved out of the macro.
 *
 * \param[in, out] xa In figure 3 of the link, would be either X6 or X5 or X1.
 *
 * \param[in, out] xb In figure 3 of the link, would be either X2 or X3 or X7.
 *
 * \param[in, out] k1 In figure 3 of the link, would be either alpha or delta or
 * êta.
 *
 * \param[in, out] k2 In figure 3 of the link, would be either bêta or
 * epsilon or thêta.
 *
 * \param[in, out] sh, the shift operator, since this is fixed
 * point arithmetic, the results might need to be shift back to the correct
 * value to continue the calculations.
 */
#define xmul(xa, xb, k1, k2, sh)                               \
  n = k1 * (xa + xb), p = xa, xa = (n + (k2 - k1) * xb) >> sh, \
  xb = (n - (k2 + k1) * p) >> sh

/**
 * \fn unsigned char NextBit(unsigned char *file_content, unsigned int &index,
 * unsigned char &bit_index)
 *
 * \brief Returns the next bit in the stream.
 *
 * This function extract the bit at the specified position in the stream. When
 * the bit_index value reach 0, the bit_index is set to 8 and the index is
 * incremented by 1. If a value outside of [1 ; 8] is provided for the
 * bit_index, the function will throw an out_of_range error.
 *
 * If a 0xFF byte if found in the stream, then, if the following byte is not
 * 0x00, the bit is extracted, else an error is thrown if the DNL byte is not
 * found. If 0xFF00 is encountered, when accessing the lowest bit in the 0xFF,
 * the index will be incremented by 2 instead of 1.
 *
 * \param[in] file_content Unsigned char pointer pointing to the data of the
 * JPEG file.
 *
 * \param[in, out] index The position of the cursor reading the file.
 *
 * \param[in, out] bit_index The index of the bit to extract. A value of 8
 * extract the highest bit in the byte, a value of 1 extract the lowest.
 *
 * \return Either 0 or 1, the value of the bit.
 */
unsigned char NextBit(unsigned char *file_content, unsigned int &index,
                      unsigned char &bit_index) {
  unsigned char current_byte = file_content[index], bit;
  if (bit_index < 1 || bit_index > 8) {
    throw std::out_of_range(
        "The index provided for the bit is out of range with a value of " +
        bit_index);
  }
  if (bit_index == 8) {
    if (current_byte == 0xFF) {
      if (!(file_content[index + 1] == 0x00)) {
        if (file_content[index + 1] == DEFINE_NUMBER_OF_LINE) {
          throw std::runtime_error("DNL byte found");
        } else {
          throw std::runtime_error(
              "Error while parsing the file, DNL byte expected but something "
              "else found");
        }
      }
    }
  }

  bit = current_byte << (8 - bit_index);
  bit = bit >> 7;
  bit_index -= 1;

  if (bit_index == 0) {
    if (current_byte == 0xFF) {
      index += 2;
    } else {
      index += 1;
    }
    bit_index = 8;
  }
  return bit;
}

/**
 * \fn void IDCT(int *new_block)
 *
 * \brief Calculate the inverse cosinus transform over the block provided. Naive
 * implementation.
 *
 * This function is to be rewritten to fit the vector calculations.
 *
 * \param[in, out] new_block A pointer to the data to modify.
 *
 * \deprecated
 */
void IDCT(std::vector<int> new_block) {
  float result, cu, cv;
  int temp_operation[64];

  for (size_t c = 0; c < 8; c++) {
    for (size_t r = 0; r < 8; r++) {
      temp_operation[c * 8 + r] = new_block[c * 8 + r];
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
          result += cu * cv * temp_operation[v + u * 8] *
                    cos((2 * x + 1) * u * M_PI / 16.0) *
                    cos((2 * y + 1) * v * M_PI / 16.0);
        }
      }
      result = result / 4 + 128;
      if (result > 255) {
        new_block[y + x * 8] = 255;
      } else if (result < 0) {
        new_block[y + x * 8] = 0;
      } else {
        new_block[y + x * 8] = (int)result;
      }
    }
  }
}

/**
 * \fn void FastIDCT1D(std::vector<int> &x, std::vector<int> &y, const int
 x_offset, const int y_offset, const int ps, const int half, const int
 y_line_length)
 *
 * \brief Compute the Fast one dimension IDCT. This implementation was taken
 from <a href="http://halicery.com/Image/idct.html">here</a> and is the
 inverse implementation of the algorithm developped by Loeffler et al.
 *
 * Since the function receives a vector containing all the image coefficients as
 input, the function uses offsets to position the cursor at the correct place.
 * The output of this function is actually: 2 * sqrt(2) * real_coefficient, it
 is left like that because, for a 2D IDCT, the factor will be of 2 * sqrt(2) * 2
 * sqrt(2) = 8, thus only requiring a shift of three to the right to divide.
 *
 *
 * \param[in,out] x A vector containing the input data.
 *
 * \param[in,out] y A vector to hold the output data.
 *
 * \param[in] x_offset The offset for the input, should place the cursor at the
 first component in the input vector.
 *
 * \param[in] y_offset The offset for the input, should place the cursor at the
 first component in the output vector. The results will be stored in the column
 of the output.
 *
 * \param[in] ps This is the precision factor, for 1D IDCT it should be set to
 9, for a 2D IDCT, should be set to 12 (9 for the precision + 3 for the final
 shift).
 *
 * \param[in] half The position of the "comma" during the IDCT, should be 1 less
 than ps. It will be used for rounding. Set to 0 to have no rounding at all.
 *
 * \param[in] y_line_length The length of a line of the image. The data is
 stored on one array, line after line, to get a 8*8 block, this is required.
 */
void FastIDCT1D(std::vector<int> &x, std::vector<int> &y, const int x_offset,
                const int y_offset, const int ps, const int half,
                const int y_line_length) {
  int p, n;
  x.at(x_offset + 0) <<= 9, x.at(x_offset + 1) <<= 7, x.at(x_offset + 3) *= 181,
      x.at(x_offset + 4) <<= 9, x.at(x_offset + 5) *= 181, x.at(x_offset + 7) <<= 7;
  xmul(x.at(x_offset + 6), x.at(x_offset + 2), 277, 669, 0);
  xadd3(x.at(x_offset + 0), x.at(x_offset + 4), x.at(x_offset + 6), x.at(x_offset + 2));
  xadd3(x.at(x_offset + 1), x.at(x_offset + 7), x.at(x_offset + 3), x.at(x_offset + 5));
  xmul(x.at(x_offset + 5), x.at(x_offset + 3), 251, 50, 6);
  xmul(x.at(x_offset + 1), x.at(x_offset + 7), 213, 142, 6);
  // We get the final results, round up if required with the addition of half
  // and shift back to the result.
  y.at(y_offset + 0 * y_line_length) =
      (x.at(x_offset + 0) + x.at(x_offset + 1)) + half >> ps;
  y.at(y_offset + 1 * y_line_length) =
      (x.at(x_offset + 4) + x.at(x_offset + 5)) + half >> ps;
  y.at(y_offset + 2 * y_line_length) =
      (x.at(x_offset + 2) + x.at(x_offset + 3)) + half >> ps;
  y.at(y_offset + 3 * y_line_length) =
      (x.at(x_offset + 6) + x.at(x_offset + 7)) + half >> ps;
  y.at(y_offset + 4 * y_line_length) =
      (x.at(x_offset + 6) - x.at(x_offset + 7)) + half >> ps;
  y.at(y_offset + 5 * y_line_length) =
      (x.at(x_offset + 2) - x.at(x_offset + 3)) + half >> ps;
  y.at(y_offset + 6 * y_line_length) =
      (x.at(x_offset + 4) - x.at(x_offset + 5)) + half >> ps;
  y.at(y_offset + 7 * y_line_length) =
      (x.at(x_offset + 0) - x.at(x_offset + 1)) + half >> ps;
}

/**
 * \fn void FastIDCT2D(std::vector<int> &image, const int start_line, const int
 * start_column, const int line_length)
 *
 * \brief Compute the a fast IDCT in dimension two. The algorithm is based upon
 * the DCT implementation of Loeffler. And was taken from <a
 * href="http://halicery.com/Image/idct.html">here</a>
 *
 * For now, the function shift back to signed values and clip the results if
 * required, this is to be moved elsewere.
 *
 * \param[in, out] image A vector containing the coefficients of one component
 of the image.
 *
 * \param[in] start_line The poisition (line) of the block being processed in
 the image.
 *
 * \param[in] start_column The poisition (column) of the block being processed
 in the image.
 *
 * \param[in] line_length The length of a line in the image.
 *
 */
void FastIDCT2D(std::vector<int> &image, const int start_line,
                const int start_column, const int line_length) {
  int i;
  std::vector<int> b2(64);

  for (i = 0; i < 8; i++)
    FastIDCT1D(image, b2,
               start_line * line_length + start_column + i * line_length, i, 9,
               1 << 8, 8);  // row
  for (i = 0; i < 8; i++)
    FastIDCT1D(b2, image, i * 8, start_line * line_length + start_column + i,
               12, 1 << 11, line_length);  // col
}

void IDCTRow(std::vector<int> &x, const int x_offset, const int y_offset,
              const int line_length) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  if (!((x1 = x.at(x_offset + 4) << 11) | (x2 = x.at(x_offset + 6)) |
        (x3 = x.at(x_offset + 2)) | (x4 = x.at(x_offset + 1)) |
        (x5 = x.at(x_offset + 7)) | (x6 = x.at(x_offset + 5)) |
        (x7 = x.at(x_offset + 3)))) {
    x.at(x_offset + 0) = x.at(x_offset + 1) = x.at(x_offset + 2) = x.at(x_offset + 3) =
        x.at(x_offset + 4) = x.at(x_offset + 5) = x.at(x_offset + 6) = x.at(x_offset + 7) =
            x.at(x_offset + 0) << 3;
    return;
  }

  x0 = (x.at(x_offset + 0) << 11) +
       128; /* for proper rounding in the fourth stage */

  /* first stage */
  x8 = W7 * (x4 + x5);
  x4 = x8 + (W1 - W7) * x4;
  x5 = x8 - (W1 + W7) * x5;
  x8 = W3 * (x6 + x7);
  x6 = x8 - (W3 - W5) * x6;
  x7 = x8 - (W3 + W5) * x7;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2);
  x2 = x1 - (W2 + W6) * x2;
  x3 = x1 + (W2 - W6) * x3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  x.at(x_offset + 0) = (x7 + x1) >> 8;
  x.at(x_offset + 1) = (x3 + x2) >> 8;
  x.at(x_offset + 2) = (x0 + x4) >> 8;
  x.at(x_offset + 3) = (x8 + x6) >> 8;
  x.at(x_offset + 4) = (x8 - x6) >> 8;
  x.at(x_offset + 5) = (x0 - x4) >> 8;
  x.at(x_offset + 6) = (x3 - x2) >> 8;
  x.at(x_offset + 7) = (x7 - x1) >> 8;
}

void IDCTCol(std::vector<int> &x, const int x_offset, const int y_offset,
              const int line_length) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;
  int iclip[1024]; /* clipping table */
  int *iclp;
  iclp = iclip + 512;
  for (int i = -512; i < 512; i++)
    iclp[i] = (i < -256) ? -256 : ((i > 255) ? 255 : i);

  if (!((x1 = (x.at(x_offset + line_length * 4) << 8)) |
        (x2 = x.at(x_offset + line_length * 6)) |
        (x3 = x.at(x_offset + line_length * 2)) |
        (x4 = x.at(x_offset + line_length * 1)) |
        (x5 = x.at(x_offset + line_length * 7)) |
        (x6 = x.at(x_offset + line_length * 5)) |
        (x7 = x.at(x_offset + line_length * 3)))) {
    x.at(x_offset + line_length * 0) = x.at(x_offset + line_length * 1) =
        x.at(x_offset + line_length * 2) = x.at(x_offset + line_length * 3) =
            x.at(x_offset + line_length * 4) = x.at(x_offset + line_length * 5) =
                x.at(x_offset + line_length * 6) = x.at(x_offset + line_length * 7) =
                    iclp[(x.at(x_offset + line_length * 0) + 32) >> 6];
    return;
  }

  x0 = (x.at(x_offset + line_length * 0) << 8) + 8192;

  /* first stage */
  x8 = W7 * (x4 + x5) + 4;
  x4 = (x8 + (W1 - W7) * x4) >> 3;
  x5 = (x8 - (W1 + W7) * x5) >> 3;
  x8 = W3 * (x6 + x7) + 4;
  x6 = (x8 - (W3 - W5) * x6) >> 3;
  x7 = (x8 - (W3 + W5) * x7) >> 3;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2) + 4;
  x2 = (x1 - (W2 + W6) * x2) >> 3;
  x3 = (x1 + (W2 - W6) * x3) >> 3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  x.at(x_offset + line_length * 0) = iclp[(x7 + x1) >> 14];
  x.at(x_offset + line_length * 1) = iclp[(x3 + x2) >> 14];
  x.at(x_offset + line_length * 2) = iclp[(x0 + x4) >> 14];
  x.at(x_offset + line_length * 3) = iclp[(x8 + x6) >> 14];
  x.at(x_offset + line_length * 4) = iclp[(x8 - x6) >> 14];
  x.at(x_offset + line_length * 5) = iclp[(x0 - x4) >> 14];
  x.at(x_offset + line_length * 6) = iclp[(x3 - x2) >> 14];
  x.at(x_offset + line_length * 7) = iclp[(x7 - x1) >> 14];
}

/* two dimensional inverse discrete cosine transform */
void FastIDCT2D_Second(std::vector<int> &image, const int start_line,
                       const int start_column, const int line_length) {

  for (int i = 0; i < 8; i++)
    IDCTRow(image, start_line * line_length + start_column + i * line_length,
             i, 8);

  for (int i = 0; i < 8; i++)
    IDCTCol(image, start_line * line_length + start_column + i, 0,
             line_length);
}

void DeLevelShift(JPEGImage &image) {
  for (int component = 0; component < image.GetNumberOfComponent();
       component++) {
    std::pair<int, int> size = image.GetComponentShape(component);
    for (size_t row = 0; row < size.first; row++) {
      for (size_t col = 0; col < size.second; col++) {
        image.at(row, col, component) += 128;
        if (image.at(row, col, component) > 255) {
          image.at(row, col, component) = 255;
        } else if (image.at(row, col, component) < 0) {
          image.at(row, col, component) = 0;
        }
      }
    }
  }
}

/**
 * \fn void YCbCrToBGR(JPEGImage &image, std::vector<int> shape)
 *
 * \brief This function converts a JPEGImage from YCbCr To BGR.
 *
 * \param[in,out] image The JPEGImage to have its space changed.
 *
 * \param[in] shape The shape of the image (argument to be removed).
 *
 */
void YCbCrToBGR(JPEGImage &image, std::vector<int> shape) {
  int R, G, B;
  double alpha;
  int rows = shape[0], cols = shape[1];
  int num_component = shape[2];
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < cols; col++) {
      
      R = image.at(row, col, 0) + 1.402 * (image.at(row, col, 2) - 128);
      G = image.at(row, col, 0) - 0.34414 * (image.at(row, col, 1) - 128) -
          0.71414 * (image.at(row, col, 2) - 128);
      B = image.at(row, col, 0) + 1.772 * (image.at(row, col, 1) - 128);
      if (num_component == 4) {
        alpha = image.at(row, col, 3) / 255.0;
      }

      if (R > 255) {
        R = 255;
      }
      if (G > 255) {
        G = 255;
      }
      if (B > 255) {
        B = 255;
      }
      if (R < 0) {
        R = 0;
      }
      if (G < 0) {
        G = 0;
      }
      if (B < 0) {
        B = 0;
      }

      if (num_component == 4) {
        // Wat ? Yolo, formula not in the norm.
        image.at(row, col, 2) = alpha * (255 - R);
        image.at(row, col, 1) = alpha * (255 - G);
        image.at(row, col, 0) = alpha * (255 - B);
      } else {
        image.at(row, col, 2) = R;
        image.at(row, col, 1) = G;
        image.at(row, col, 0) = B;
      }
    }
  }
}