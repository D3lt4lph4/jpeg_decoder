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
  // Decoder related functions.
  void InitializeDecoder();
  void DecoderSetup();

  // Decoding functions
  void DecodeFrame(unsigned char encoding_process_type);
  void DecodeScan(unsigned char encoding_process_type);
  void DecodeRestartIntervalProgressive();
  void DecodeMCUProgressive();
  void DecodeRestartIntervalBaseline();

  // Parsing function
  void ParseQuantizationTable();
  void InterpretFrameHeader(unsigned char encoding_process_type);
  void ParseHuffmanTableSpecification();
  void ParseScanHeader();
  void ParseComment();
  void ParseJFIFSegment();
  void ProcessData();

  // Baseline functions
  void ResetDecoderBaseline();
  unsigned char DecodeBaseline(HuffmanTable used_table);
  int ReceiveBaseline(unsigned char decoded_dc);
  int ExtendedBaseline(int diff, unsigned char ssss);

  // Progressive encoding functions
  void ResetDecoderProgressive();

  // Utility function required in the norm.
  unsigned char NextBit();
  void DecodeACCoefficients(cv::Mat *new_block, HuffmanTable used_table);
  unsigned char DecodeZZ(unsigned char ssss);

  // Functions to process the huffman table from the image stream.
  std::vector<unsigned char> GenerateSizeTable(std::vector<unsigned char> bits);
  std::vector<unsigned short> GenerateCodeTable(std::vector<unsigned char> huffsize);
  void DecoderTables(HuffmanTable *table_processed);

  // Marker related functions
  bool IsMarker();
  unsigned char *GetMarker();

  unsigned int restart_interval;
  int current_define_quantization_table_, data_unit_per_mcu_, decoding_level_;
  cv::Mat current_image_, current_thumbnail_;
  std::string current_filename_;

  // The different tables
  std::map<unsigned char, QuantizationTable> quantization_tables_;
  std::map<unsigned char, HuffmanTable> dc_huffman_tables_, ac_huffman_tables_;

  // The headers
  FrameHeader current_frame_header_;
  ScanHeader current_scan_;
  JFIFHeader current_jfif_header;

  // File crossing variable.
  unsigned char *current_file_content_, bit_index;
  int current_index_;

  // Variables required for the functions in the norm
  char next_bit_count_;
  char last_k_;
};

#endif