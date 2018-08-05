#ifndef __JPEG_HUFFMAN_DECODER__
#define __JPEG_HUFFMAN_DECODER__

#include <opencv2/core.hpp>

#include "JPEGType.hpp"

// Functions to parse the huffman table as found in the stream.
std::pair<unsigned char, std::vector<unsigned char>> GenerateSizeTable(
    std::vector<unsigned char>);
std::vector<unsigned char> GenerateCodeTable(std::vector<unsigned char>);
std::tuple<std::vector<unsigned char>, std::vector<unsigned char>,
           std::vector<unsigned char>>
DecoderTables(HuffmanTable *table_being_parsed);

// Functions to process the stream
unsigned char Decode(unsigned char *stream, unsigned int *index,
                     unsigned char *bit_index, HuffmanTable used_table);
int Receive(unsigned char number_of_bits, unsigned char *stream,
            unsigned int *index, unsigned char *bit_index);
int Extended(int diff, unsigned char number_of_bits);
void DecodeACCoefficients(unsigned char *stream, unsigned int *index,
                          unsigned char *bit_index, cv::Mat *new_block,
                          HuffmanTable used_table);
unsigned char DecodeZZ(unsigned char *stream, unsigned int *index,
                       unsigned char *bit_index, unsigned char ssss);

#endif