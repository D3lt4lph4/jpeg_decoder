#ifndef __JPEG_DECODER__
#define __JPEG_DECODER__

#include <map>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGType.hpp"

class JPEGDecoder
{
public:
  cv::Mat Decode(std::string file_name, int level);
  friend std::ostream &operator<<(std::ostream &out,
                                  const JPEGDecoder &decoder);

private:
  void DecodeHuffman();
  void DecodeRunLengthEncoding();
  void InverseQuantification();
  void InverseDirectCosineTransform();

  unsigned char *GetMarker(unsigned char *file_content, int *index,
                           int size = 1, bool FF_expected = true);
  void ParseQuantizationTable(unsigned char *file_content, int *index);
  void ParseFrameHeader(unsigned char *file_content, int *index, unsigned char encoding_process_type);
  void ParseHuffmanTableSpecification(unsigned char *file_content, int *index);
  void ParseScanHeader(unsigned char *file_content, int *index);
  void ParseComment(unsigned char *file_content, int *index);
  bool GetFileInformation(unsigned char *file_content, int *index);
  void ProcessData(unsigned char *file_content, int *index);

  // Function required in the norm
  char NextBit(unsigned char *file_content, int *index);
  char Decode(unsigned char *file_content, int *index, HuffmanTable used_table);

  std::vector<char> GenerateSizeTable(std::vector<char> bits);
  std::vector<int> GenerateCodeTable(std::vector<char> huffsize);

  int current_version_, current_unit_, horizontal_pixel_density_,
      vertical_pixel_density_, thumbnail_horizontal_pixel_count_,
      thumbnail_vertical_pixel_count_, current_define_quantization_table_;
  cv::Mat current_image_, current_thumbnail_;
  std::string current_filename_;

  std::map<unsigned char, QuantizationTable> quantization_tables_;
  std::map<unsigned char, HuffmanTable> huffman_tables_;

  FrameHeader current_frame_header_;
  ScanHeader current_scan_;

  // Variables required for the functions in the norm
  char next_bit_count_;
  char last_k_;
};

#endif