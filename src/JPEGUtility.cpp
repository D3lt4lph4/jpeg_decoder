#include <stdexcept>

#include "JPEGType.hpp"
#include "JPEGUtility.hpp"

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
