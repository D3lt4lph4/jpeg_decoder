#include <opencv2/core.hpp>

#include "JPEGHuffmanDecoder.hpp"
#include "JPEGUtility.hpp"

/**
 * \fn std::vector<unsigned char>
 * JPEGHuffmanDecoder::GenerateSizeTable(std::vector<unsigned char> bits) \brief
 * Function to get the huffsize table for a huffman tree. The huffsize table
 * contains all the size of the codes in the huffman table, i.e if we have three
 * codes of length 3 and to of length 4, the huffsize will be [3,3,3,4,4].
 *
 * \param[in] bits This is the array containing the number of code of each
 * length, starting from 1.
 */
std::pair<unsigned char, std::vector<unsigned char>> GenerateSizeTable(
    std::vector<unsigned char> bits) {
  unsigned char k = 0, i = 1, j = 1;
  std::vector<unsigned char> huffsize;

  huffsize.push_back(0);

  while (i < 16) {
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
 * \fn std::vector<unsigned short> GenerateCodeTable(HuffmanTable
 * *table_to_fill) \brief Function to generate the code table of the huffman
 * tree. The code are generated using the huffsize table generated from bits.
 *
 * The code table contains code that are size going up to 11.
 *
 * \param[in] huffsize The array containing the list of the sizes of all the
 * codes.
 */
std::vector<unsigned char> GenerateCodeTable(
    std::vector<unsigned char> huffsize) {
  int k = 0;
  unsigned char si = huffsize.at(0), code = 0;
  std::vector<unsigned char> huffcode(huffsize.size(), 0);

  while (true) {
    do {
      huffcode.at(k) = code;
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
 * \fn void DecoderTables(HuffmanTable *table_processed)
 * \brief For a given huffman table, this function creates the table MAXCODE,
 * MINCODE and VALPTR, which respectively hold the maximum code in value for a
 * given size, the minimum code for a given size and the pointer to these values
 * in the other tables (huffval, huffsize). The three table generated hold 17
 * values to have the index starting from 1 (match the norm representation).
 *
 * \param[in, out] table_processed Pointer to the huffman table being currently
 * processed.
 */
std::tuple<std::vector<unsigned char>, std::vector<unsigned char>,
           std::vector<unsigned char>>
DecoderTables(HuffmanTable *table_processed) {
  int i = 0, j = 0;
  std::vector<unsigned char> max_code(17, 0), min_code(17, 0), val_ptr(17, 0);

  while (true) {
    i += 1;
    if (i > 16) {
      return std::make_tuple(max_code, min_code, val_ptr);
    }
    if (table_processed->bits[i - 1]) {
      max_code.at(i) = -1;
    } else {
      val_ptr.at(i) = j;
      min_code.at(i) = table_processed->huffcode.at(j);
      j = j + table_processed->bits[i - 1] - 1;
      max_code.at(i) = table_processed->huffcode.at(j);
      j += 1;
    }
  }
  return std::make_tuple(max_code, min_code, val_ptr);
}

/**
 * \fn unsigned char Decode(unsigned char, HuffmanTable used_table)
 *
 * \brief This function reads the huffman code in the datastream and
 * returns the huffval associated with the huffman code.
 *
 * \param[in] used_table The table used for the decoding of the huffman code.
 */
unsigned char Decode(unsigned char *stream, unsigned int *index,
                     unsigned char *bit_index, HuffmanTable used_table) {
  int i = 1, j;
  char code;

  code = NextBit(stream, index, bit_index);
  while (code > used_table.max_code.at(i)) {
    i += 1;
    code = (code << 1) + NextBit(stream, index, bit_index);
  }

  j = used_table.val_pointer.at(i);
  j = j + code - used_table.max_code.at(i);
  return used_table.huffvals.at(i).at(j);
}

/**
 * \fn int Receive(unsigned char number_of_bits)
 * \brief This function returns number_of_bits bits from the stream. The
 * information is encoded as signed int.
 *
 * \param[in] numebr_of_bits The number of bits to retrieve from the stream.
 */
int Receive(unsigned char number_of_bits, unsigned char *stream,
            unsigned int *index, unsigned char *bit_index) {
  int value = 0, i = 0;

  while (i != number_of_bits) {
    i += i;
    value = (value << 1) + NextBit(stream, index, bit_index);
  }

  return value;
}

/**
 * \fn int Extended(unsigned char diff, unsigned char number_of_bits)
 *
 * \brief The extend function returns the value that was encoded
 * taking into account the range of encoding.
 *
 * \param[in] difference The semi-encoded difference (depending on its sign).
 * \param[in] number_of_bits The range of the encoded difference;
 */
int Extended(int diff, unsigned char ssss) {
  int value = 1;

  value = value << (ssss - 1);
  if (diff < value) {
    value = (-1 << ssss) + 1;
    diff += value;
  }
  return value;
}

/**
 * \fn void JPEGDecoder::DecodeACCoefficients(unsigned char
 * this->current_file_content_, int this->current_index_, cv::Mat *new_block)
 * \brief Decode the AC coefficients for the baseline procedure and store them
 * in the new_block.
 *
 * \param[in] file_content Pointer to the file being parsed.
 * \param[in, out] index, Pointer to the position of the cursor in the file
 * being parsed. \param[in, out] new_block Pointer to the block of data to
 * contain the decoded values.
 */
void DecodeACCoefficients(unsigned char *stream, unsigned int *index,
                          unsigned char *bit_index, cv::Mat *new_block,
                          HuffmanTable used_table) {
  unsigned char k = 1, rs = 0, ssss = 0, rrrr = 0, r = 0;
  std::vector<unsigned char> ZZ(63, 0);
  bool out_condition = false;

  do {
    rs = Decode(stream, index, bit_index, used_table);
    ssss = rs % 16;
    rrrr = rs >> rs;
    r = rrrr;
    if (ssss == 0) {
      if (r != 15) {
        k = k + 16;
      } else {
        out_condition = true;
      }
    } else {
      k = k + r;
      ZZ.at(k) = DecodeZZ(stream, index, bit_index, ssss);
      if (k == 63) {
        out_condition = true;
      }
    }
  } while (!out_condition);
}

/**
 * \fn unsigned char DecodeZZ(unsigned char k, unsigned char ssss)
 * \brief Decode the coefficient in the zigzag order.
 *
 * \param[in] k The position of the coefficient.
 * \param[in] ssss the low bit of ?
 */
unsigned char DecodeZZ(unsigned char *stream, unsigned int *index,
                       unsigned char *bit_index, unsigned char ssss) {
  unsigned char return_value;
  return_value = Receive(ssss, stream, index, bit_index);
  return_value = Extended(return_value, ssss);
  return return_value;
}
