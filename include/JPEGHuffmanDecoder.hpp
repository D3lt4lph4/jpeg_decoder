#ifndef __JPEG_HUFFMAN_DECODER__
#define __JPEG_HUFFMAN_DECODER__

#include <opencv2/core.hpp>

#include "JPEGType.hpp"

// Functions to parse the huffman table as found in the stream.
std::pair<unsigned char, std::vector<unsigned char>> GenerateSizeTable(
    std::vector<unsigned char>);
std::vector<unsigned short> GenerateCodeTable(std::vector<unsigned char>);
std::tuple<std::vector<int>, std::vector<int>, std::vector<unsigned char>>
DecoderTables(std::vector<unsigned char> bits,
              std::vector<unsigned short> huffcode);

// Functions to process the stream
unsigned char Decode(unsigned char *stream, unsigned int *index,
                     unsigned char *bit_index, HuffmanTable used_table);
int Receive(unsigned char number_of_bits, unsigned char *stream,
            unsigned int *index, unsigned char *bit_index);
int Extended(int diff, unsigned char number_of_bits);
std::vector<int> DecodeACCoefficients(unsigned char *stream,
                                      unsigned int *index,
                                      unsigned char *bit_index,
                                      HuffmanTable used_table);
unsigned char DecodeZZ(unsigned char *stream, unsigned int *index,
                       unsigned char *bit_index, unsigned char ssss);

#endif