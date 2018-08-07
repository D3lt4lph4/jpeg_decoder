#ifndef __JPEG_TYPE__
#define __JPEG_TYPE__

#include <map>
#include <vector>

const unsigned char FRAME_TYPE_BASELINE_DTC = 0xB0;
const unsigned char FRAME_TYPE_PROGRESSIVE = 0xB2;

const unsigned char START_OF_FRAME_BASELINE = 0xC0;
const unsigned char START_OF_FRAME_PROGRESSIVE = 0xC2;
const unsigned char DEFINE_HUFFMAN_TABLE = 0xC4;

const unsigned char START_OF_IMAGE = 0xD8;
const unsigned char END_OF_IMAGE = 0xD9;
const unsigned char START_OF_SCAN = 0xDA;
const unsigned char DEFINE_QUANTIZATION_TABLE = 0xDB;
const unsigned char DEFINE_NUMBER_OF_LINE = 0xDC;
const unsigned char DEFINE_RESTART_INTERVAL = 0xDD;
const unsigned char DEFINE_HIERARCHICAL_PROGRESSION = 0xDE;
const unsigned char EXPAND_REFERENCE_COMPONENTS = 0xDF;

const unsigned char COMMENT = 0xFE;

const unsigned char APPO = 0xE0;

const unsigned char JFIF[] = {0x4a, 0x46, 0x49, 0x46, 0x00};

struct QuantizationTable {
  unsigned char pq_;
  std::vector<unsigned short> qks_;
};

struct FrameComponentSignification {
  unsigned char horizontal_sampling_factor_, vertical_sampling_factor,
      quantization_table_selector;
};

struct JFIFHeader {
  unsigned short current_version_;
  unsigned char current_unit_, horizontal_pixel_density_,
      vertical_pixel_density_, thumbnail_horizontal_pixel_count_,
      thumbnail_vertical_pixel_count_;
};

struct FrameHeader {
  unsigned char encoding_process_type_, sample_precision_;
  unsigned int number_of_lines_, number_of_samples_per_line_,
      number_of_image_component;
  std::map<unsigned char, std::vector<unsigned char>>
      component_signification_parameters_;
};

struct ScanHeader {
  unsigned char number_of_image_components_, start_of_spectral_selection_,
      end_of_spectral_selection_, approximation_high_bit_,
      approximation_low_bit_;
  std::map<unsigned char, std::pair<unsigned char, unsigned char>>
      scan_components_specification_parameters_;
};

struct HuffmanTable {
  unsigned char table_class_, last_k_;
  std::vector<unsigned char> bits, huffsize, val_pointer, huffvals;
  std::vector<unsigned short> huffcode;
  std::vector<int> min_code, max_code;
};
#endif