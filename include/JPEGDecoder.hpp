/**
 * \file JPEGDecoder.hpp
 */

#ifndef __JPEG_DECODER__
#define __JPEG_DECODER__

#include <map>

#include "JPEGHuffmanDecoder.hpp"
#include "JPEGParser.hpp"
#include "JPEGType.hpp"

/**
 * \class JPEGDecoder
 */
class JPEGDecoder {
 public:
  // Constructors & Destructor
  JPEGDecoder();
  JPEGDecoder(unsigned char logging_level_);
  ~JPEGDecoder();

  // Class functions
  void *DecodeFile(std::string file_name, unsigned int *image_size_x,
                   unsigned int *image_size_y, int level);
  friend std::ostream &operator<<(std::ostream &out,
                                  const JPEGDecoder &decoder);

 private:
  // Decoder related functions.
  void InitializeDecoder();
  void ClearVariables();
  void DecoderSetup();

  // Decoding functions
  void DecodeToLevel();
  void DecodeFrame(unsigned char encoding_process_type);
  void DecodeScan(unsigned char encoding_process_type);
  void Dequantize(int *new_block, QuantizationTable table,
                  unsigned int component_number);

  // Baseline functions
  void ResetDecoderBaseline();
  void DecodeRestartIntervalBaseline();
  void DecodeMCUBaseline(unsigned int mcu_number, unsigned int h_max,
                         unsigned int v_max, unsigned char *bit_index,
                         int *prev);

  // Marker related functions
  bool IsMarker();
  unsigned char *GetMarker();

  // Logging functions
  void InitializeLogger();

  unsigned int restart_interval;
  int current_define_quantization_table_, data_unit_per_mcu_, decoding_level_;
  void *current_image_, *current_thumbnail_;
  std::string current_filename_;

  // The different tables
  std::map<unsigned char, QuantizationTable> quantization_tables_;
  std::map<unsigned char, HuffmanTable> dc_huffman_tables_, ac_huffman_tables_;

  // The headers
  FrameHeader frame_header_;
  ScanHeader scan_header_;
  JFIFHeader current_jfif_header;

  // File crossing variable.
  unsigned char *current_file_content_, *bit_index, *last_k_;
  unsigned int *current_index_;

  int logging_level_;

  int number_of_blocks_per_line, number_of_blocks_per_column, block_index;
};

#endif