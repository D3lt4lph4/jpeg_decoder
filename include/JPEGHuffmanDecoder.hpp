#ifndef __JPEG_HUFFMAN_DECODER__
#define __JPEG_HUFFMAN_DECODER__

#include "JPEGType.hpp"

// Functions to parse the huffman table as found in the stream.
std::pair<unsigned char, std::vector<unsigned char>> GenerateSizeTable(
    const std::vector<unsigned char> &bits);
std::vector<unsigned short> GenerateCodeTable(const std::vector<unsigned char> &huffsize);
std::tuple<std::vector<int>, std::vector<int>, std::vector<unsigned char>>
DecoderTables(const std::vector<unsigned char> &bits,
              const std::vector<unsigned short> &huffcode);

// Functions to process the stream
unsigned char Decode(unsigned char *stream, unsigned int &index,
                     unsigned char &bit_index, const HuffmanTable &used_table);
int Receive(const unsigned char number_of_bits, unsigned char *stream,
            unsigned int &index, unsigned char &bit_index);
int Extended(int diff, const unsigned char number_of_bits);
std::vector<int> DecodeACCoefficients(unsigned char *stream,
                                      unsigned int &index,
                                      unsigned char &bit_index,
                                      const HuffmanTable &used_table);
int DecodeZZ(unsigned char *stream, unsigned int &index,
             unsigned char &bit_index, const unsigned char ssss);

#endif