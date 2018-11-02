/**
 * \file JPEGHuffmanDecoder.cpp
 * \author Benjamin Deguerre
 * \brief Contains all the functions
 */

#include <iostream>

#include "JPEGHuffmanDecoder.hpp"
#include "JPEGUtility.hpp"

/**
 * \fn std::pair<unsigned char, std::vector<unsigned char>>
 * GenerateSizeTable(const std::vector<unsigned char> &bits)
 *
 * \brief Function to get the huffsize table for a huffman tree. The huffsize
 * table contains all the size of the codes in the huffman table, i.e if we have
 * three codes of length 3 and to of length 4, the huffsize will be
 * [3,3,3,4,4,0]. The zero at the end is used as marker for other functions.
 *
 * \param[in] bits This is the array containing the number of code of each
 * length, starting from 1.
 *
 * \return std::pair<unsigned char, std::vector<unsigned char>> The first
 * element of the pair is the index of the 0 value in the array and the second
 * element is the array of length, the length of the array is sum(bits) + 1
 */
std::pair<unsigned char, std::vector<unsigned char>> GenerateSizeTable(
    const std::vector<unsigned char> &bits) {
  unsigned char k = 0, i = 1, j = 1;
  std::vector<unsigned char> huffsize;

  while (i <= 16) {
    if (j > bits.at(i - 1)) {
      i += 1;
      j = 1;
    } else {
      huffsize.push_back(i);
      k += 1;
      j += 1;
    }
  }
  huffsize.push_back(0);

  return std::make_pair(k, huffsize);
}

/**
 * \fn std::vector<unsigned short> GenerateCodeTable(const std::vector<unsigned
 * char> &huffsize)
 *
 * \brief Function to generate the code table of the huffman tree. The code are
 * generated using the huffsize table generated from bits.
 *
 * The code table contains code that have size going up to 11 for baseline (8
 * bits input precision) and 15 if the input precision for for the input is 12
 * bits.
 *
 * \param[in] huffsize The array containing the list of the sizes of all the
 * codes.
 *
 * \return A vector of shorts containing sum(huffsize) codes.
 */
std::vector<unsigned short> GenerateCodeTable(
    const std::vector<unsigned char> &huffsize) {
  int k = 0;
  unsigned char si = huffsize.at(0);
  unsigned short code = 0;
  std::vector<unsigned short> huffcode;

  while (true) {
    do {
      huffcode.push_back(code);
      code += 1;
      k += 1;
    } while (huffsize.at(k) == si);

    if (huffsize.at(k) == 0) {
      return huffcode;
    }
    do {
      code = code << 1;
      si = si + 1;
    } while (huffsize.at(k) != si);
  }
}

/**
 * \fn std::tuple<std::vector<int>, std::vector<int>, std::vector<unsigned
 * char>> DecoderTables(const std::vector<unsigned char> &bits, const
 * std::vector<unsigned short> &huffcode)
 *
 * \brief For a given huffman table, this function creates the table MAXCODE,
 * MINCODE and VALPTR, which respectively hold the maximum code in value for a
 * given size, the minimum code for a given size and the pointer to these values
 * in the other tables (huffval, huffsize).
 * The three table generated hold 17 values to have the index starting from 1
 * (match the norm representation).
 *
 * \param[in] bits The array containing the number of code of each length.
 *
 * \param[in] huffcode The array containing all the codes.
 *
 * \return A tuple containing in that order, max_code, min_code, val_ptr;
 */
std::tuple<std::vector<int>, std::vector<int>, std::vector<unsigned char>>
DecoderTables(const std::vector<unsigned char> &bits,
              const std::vector<unsigned short> &huffcode) {
  int i = 0, j = 0;
  std::vector<int> max_code(17, 0), min_code(17, 0);
  std::vector<unsigned char> val_ptr(17, 0);

  while (true) {
    i += 1;
    if (i > 16) {
      return std::make_tuple(max_code, min_code, val_ptr);
    }
    if (bits[i - 1] == 0) {
      max_code.at(i) = -1;
    } else {
      val_ptr.at(i) = j;
      min_code.at(i) = huffcode.at(j);
      j = j + bits[i - 1] - 1;
      max_code.at(i) = huffcode.at(j);
      j += 1;
    }
  }
}

/**
 * \fn unsigned char Decode(unsigned char *stream, unsigned int &index, unsigned
 * char &bit_index, const HuffmanTable &used_table)
 *
 * \brief This function reads the huffman code in the datastream and
 * returns the huffval associated with the huffman code, i.e the range of the
 * encoded value.
 *
 * \param[in] stream The pointer to the file being parsed.
 *
 * \param[in, out] index The offset of the pointer reading the file.
 *
 * \param[in, out] bit_index The offset of the bit being read in the byte.
 *
 * \param[in] used_table The table used for the decoding of the huffman code.
 *
 * \return The value associated with the huffman code.
 */
unsigned char Decode(unsigned char *stream, unsigned int &index,
                     unsigned char &bit_index, const HuffmanTable &used_table) {
  int i = 1;
  unsigned short j;
  unsigned short code;

  code = NextBit(stream, index, bit_index);
  while (code > used_table.max_code.at(i)) {
    i += 1;
    code = (code << 1) + NextBit(stream, index, bit_index);
  }

  j = used_table.val_pointer.at(i);
  j = j + code - used_table.min_code.at(i);
  return used_table.huffvals.at(j);
}

/**
 * \fn int Receive(const unsigned char number_of_bits, unsigned char *stream,
 * unsigned int &index, unsigned char &bit_index)
 *
 * \brief This function returns number_of_bits bits from the stream. The
 * information is encoded as signed int.
 *
 * \param[in] number_of_bits The number of bits to retrieve from the stream.
 *
 * \param[in] stream The pointer to the file being parsed.
 *
 * \param[in] index The offset of the pointer reading the file.
 *
 * \param[in] bit_index The offset of the bit being read in the byte.
 *
 * \return The value read in the stream.
 */
int Receive(const unsigned char number_of_bits, unsigned char *stream,
            unsigned int &index, unsigned char &bit_index) {
  int value = 0;
  unsigned char i = 0;
  unsigned char temp;

  while (i != number_of_bits) {
    i += 1;
    value = (value << 1) + NextBit(stream, index, bit_index);
  }

  return value;
}

/**
 * \fn int Extended(int diff, const unsigned char number_of_bits)
 *
 * \brief The extend function returns the value that was encoded
 * taking into account the range of encoding.
 *
 * \param[in] difference The semi-encoded difference (depending on its sign).
 * \param[in] number_of_bits The range of the encoded difference;
 *
 * \return The corrected diff value.
 */
int Extended(int diff, const unsigned char ssss) {
  int value = 1;

  value = value << (ssss - 1);
  if (diff < value) {
    value = (-1 << ssss) + 1;
    diff += value;
  }
  return diff;
}

/**
 * \fn std::vector<int> DecodeACCoefficients(unsigned char *stream, unsigned int
 * &index, unsigned char &bit_index, const HuffmanTable &used_table)
 *
 * \brief Decode the AC coefficients for the baseline procedure and store them
 * in the new_block.
 *
 * \param[in] stream Pointer to the file being parsed.
 *
 * \param[in, out] index The position of the cursor in the file being parsed.
 *
 * \param[in, out] bit_index The offset of the bit being read in the byte.
 * 
 * \param[in] used_table The table used to decode the coefficients.
 * 
 * \return A vector containing all the AC coefficients.
 */
std::vector<int> DecodeACCoefficients(unsigned char *stream,
                                      unsigned int &index,
                                      unsigned char &bit_index,
                                      const HuffmanTable &used_table) {
  unsigned char k = 1, rs = 0, ssss = 0, rrrr = 0, r = 0;
  std::vector<int> ZZ(63, 0);
  bool out_condition = false;

  do {
    rs = Decode(stream, index, bit_index, used_table);
    ssss = rs % 16;
    rrrr = (rs & 240) >> 4;
    r = rrrr;
    if (ssss == 0) {
      if (r == 15) {
        k = k + 16;
      } else {
        out_condition = true;
      }
    } else {
      k = k + r;
      ZZ.at(k - 1) = DecodeZZ(stream, index, bit_index, ssss);
      if (k == 63) {
        out_condition = true;
      }
      k += 1;
    }
  } while (!out_condition);
  return ZZ;
}

/**
 * \fn int DecodeZZ(unsigned char *stream, unsigned int &index, unsigned char &bit_index, const unsigned char ssss)
 * 
 * \brief Decode the coefficient in the zigzag order.
 *
 * \param[in] stream Pointer to the file being parsed.
 * 
 * \param[in, out] index The position of the cursor in the file being parsed.
 * 
 * \param bit_index The offset of the bit being read in the byte.

 * \param[in] ssss The low bit of the RRRRSSSS encoding for the AC coefficients,
 * ssss represent the range of the AC coefficient and the number of bits to
 * follow representing the ac coefficient.
 * 
 * \return The decoded value from the stream.
 */
int DecodeZZ(unsigned char *stream, unsigned int &index,
             unsigned char &bit_index, const unsigned char ssss) {
  int return_value;
  return_value = Receive(ssss, stream, index, bit_index);
  return_value = Extended(return_value, ssss);
  return return_value;
}
