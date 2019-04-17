#ifndef __JPEG_DECODER__
#define __JPEG_DECODER__

#include <map>
#include <memory>

#include "JPEGHuffmanDecoder.hpp"
#include "JPEGParser.hpp"
#include "JPEGType.hpp"
#include "JPEGUtils.hpp"


// Class functions
DecodeFile(const std::string file_name, const int level, JPEGImage *);
ImageType GetType(const std::string file_name);


  // Decoder related functions.
  void ClearVariables();

  // Decoding functions
  void DecodeToLevel();
void DecodeFrame(const unsigned char encoding_process_type);
ImageType DecodeFrameType(const unsigned char encoding_process_type);
void DecodeScan(const unsigned char encoding_process_type);
void Dequantize(const int component_number, const int start_row,
                const int start_col, const QuantizationTable &table);
  void Upscale();

  // Baseline functions
void DecodeRestartIntervalBaseline();
void DecodeMCUBaseline(const unsigned int mcu_number,
                        const unsigned int h_max, const unsigned int v_max,
                        unsigned char &bit_index, std::vector<int> &prev);

// Marker related functions
bool IsMarker();
std::unique_ptr<unsigned char> GetMarker();


  unsigned int restart_interval;
  int current_define_quantization_table_, data_unit_per_mcu_, decoding_level_;
  std::string current_filename_;

  // The different tables
  std::map<unsigned char, QuantizationTable> quantization_tables_;
  std::map<unsigned char, HuffmanTable> dc_huffman_tables_, ac_huffman_tables_;

  // The headers
  FrameHeader frame_header_;
  ScanHeader scan_header_;
  JFIFHeader current_jfif_header;

  // File crossing variable.
  unsigned char *current_file_content_, bit_index, last_k_;
  unsigned int current_index_;


  int number_of_blocks_per_line, number_of_blocks_per_column, block_index;

#endif