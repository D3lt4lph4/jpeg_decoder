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

  // Setup Functions
  void InitializeDecoder();
  void DecoderSetup();
  void ResetDecoderProgressive();
  void ResetDecoderBaseline();

  // Decoding functions
  void DecodeFrame(unsigned char *file_content, int *index, unsigned char encoding_process_type);
  void DecodeScan(unsigned char *file_content, int *index, unsigned char encoding_process_type);
  void DecodeRestartIntervalProgressive(unsigned char *file_content, int *index);
  void DecodeMCUProgressive(unsigned char *file_content, int *index);
  void DecodeRestartIntervalBaseline(unsigned char *file_content, int *index);
  void DecodeMCUBaseline(unsigned char *file_content, int *index);

  // Parsing function
  void ParseQuantizationTable(unsigned char *file_content, int *index);
  void InterpretFrameHeader(unsigned char *file_content, int *index, unsigned char encoding_process_type);
  void ParseHuffmanTableSpecification(unsigned char *file_content, int *index);
  void ParseScanHeader(unsigned char *file_content, int *index);
  void ParseComment(unsigned char *file_content, int *index);
  void ParseJFIFSegment(unsigned char *file_content, int *index);
  void ProcessData(unsigned char *file_content, int *index);

  // Function required in the norm
  unsigned char NextBit(unsigned char *file_content, int *index);
  unsigned char DecodeBaseline(unsigned char *file_content, int *index, HuffmanTable used_table);
  int ReceiveBaseline(unsigned char decoded_dc);
  int ExtendedBaseline(unsigned char diff, unsigned char decoded_dc);
  void DecodeACCoefficients(unsigned char *file_content, int *index);

  std::vector<unsigned char> GenerateSizeTable(std::vector<unsigned char> bits);
  void GenerateCodeTable(HuffmanTable *table_to_fill);

  // Marker related functions
  bool IsMarker(unsigned char *file_content, int index);
  unsigned char *GetMarker(unsigned char *file_content, int *index,
                           int size = 1, bool FF_expected = true);

  unsigned int restart_interval;
  int current_version_, current_unit_, horizontal_pixel_density_,
      vertical_pixel_density_, thumbnail_horizontal_pixel_count_,
      thumbnail_vertical_pixel_count_, current_define_quantization_table_, data_unit_per_mcu_;
  cv::Mat current_image_, current_thumbnail_;
  std::string current_filename_;

  // The different tables
  std::map<unsigned char, QuantizationTable> quantization_tables_;
  std::map<unsigned char, HuffmanTable> huffman_tables_;

  // The headers
  FrameHeader current_frame_header_;
  ScanHeader current_scan_;
  JFIFHeader current_jfif_header;

  // File crossing variable.
  unsigned char *current_file_content_;
  int current_index_;

  // Variables required for the functions in the norm
  char next_bit_count_;
  char last_k_;
};

#endif