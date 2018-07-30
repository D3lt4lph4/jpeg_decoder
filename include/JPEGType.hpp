#ifndef __JPEG_TYPE__
#define __JPEG_TYPE__

const unsigned char START_OF_IMAGE = 0xD8;
const unsigned char END_OF_IMAGE = 0xD9;
const unsigned char START_OF_FRAME_BASELINE = 0xC0;
const unsigned char START_OF_FRAME_PROGRESSIVE = 0xC2;
const unsigned char DEFINE_HUFFMAN_TABLE = 0xC4;
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
  std::vector<unsigned int> qks_;
};

struct FrameComponentSignification {
  unsigned char horizontal_sampling_factor_, vertical_sampling_factor,
      quantization_table_selector;
};

struct FrameHeader {
  unsigned char sample_precision_;
  unsigned int number_of_lines_, number_of_samples_per_line,
      number_image_component;
  std::vector<FrameComponentSignification> component_signification_parameters;
};

struct ScanHeader {
  unsigned char number_of_image_components_, ss_, se_, ah_, ai;
  std::vector<std::vector<unsigned int>>
      scan_component_specification_parameters_;
};


struct HuffmanTable {
  unsigned char tc_, th_;
  std::vector<unsigned int> number_of_codes_of_length_i_;
  std::vector<std::vector<unsigned int>> huffman_code_associated_values_;
};
#endif