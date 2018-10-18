#ifndef __JPEG_PARSER__
#define __JPEG_PARSER__

#include "JPEGHuffmanDecoder.hpp"
#include "JPEGType.hpp"

JFIFHeader ParseJFIFSegment(unsigned char* file_content, unsigned int* index);
std::string ParseComment(unsigned char* file_content, unsigned int* index);
FrameHeader ParseFrameHeader(unsigned char* file_content, unsigned int* index,
                             unsigned char encoding_process_type);
ScanHeader ParseScanHeader(unsigned char* file_content, unsigned int* index);
std::pair<unsigned char, QuantizationTable> ParseQuantizationTable(
    unsigned char* file_content, unsigned int* index);
std::vector<std::pair<unsigned char, HuffmanTable>>
ParseHuffmanTableSpecification(unsigned char* file_content,
                               unsigned int* index);
void ParseApplicationBlock(unsigned char* file_content, unsigned int* index);

#endif