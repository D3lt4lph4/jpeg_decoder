#include <stdexcept>

#include "JPEGType.hpp"
#include "JPEGUtility.hpp"

/**
 * \fn unsigned char NextBit()
 * \brief Returns the next bit in the stream.
 */
unsigned char NextBit(unsigned char *file_content, unsigned int *index,
                      unsigned char *bit_index) {
  unsigned char current_byte = file_content[*index], bit;

  if (*bit_index == 0) {
    *bit_index = 8;
    if (current_byte == 0xFF) {
      if (!(file_content[*index + 1] == 0x00)) {
        if (file_content[*index + 1] == DEFINE_NUMBER_OF_LINE) {
        } else {
          throw std::runtime_error(
              "Error while parsing the file, DNL byte expected but something "
              "else found");
        }
      } else {
        if ((*bit_index - 1) == 0) {
          *index += 1;
        }
      }
    }
  }

  bit = current_byte << (8 - *bit_index);
  bit = bit >> 7;
  *bit_index -= 1;

  if (*bit_index == 0) {
    *index += 1;
  }
  return bit;
}
