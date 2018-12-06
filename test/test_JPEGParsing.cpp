#include "JPEGDecoder.hpp"
#include "JPEGType.hpp"
#include "gmock/gmock.h"

TEST(ParsingTest, test_JFIF) {
  unsigned char file_content[] = {0x00, 0x10, 0x4A, 0x46, 0x49, 0x46,
                                  0x00, 0x01, 0x01, 0x01, 0x00, 0x48,
                                  0x00, 0x48, 0x00, 0x00};
  unsigned int index = 0;
  JFIFHeader header = ParseJFIFSegment(file_content, index);

  ASSERT_EQ(header.current_version_, 257);
  ASSERT_EQ(header.current_unit_, 1);
  ASSERT_EQ(header.horizontal_pixel_density_, 72);
  ASSERT_EQ(header.vertical_pixel_density_, 72);
  ASSERT_EQ(header.thumbnail_horizontal_pixel_count_, 0);
  ASSERT_EQ(header.thumbnail_vertical_pixel_count_, 0);

  ASSERT_EQ(index, 16);
}

TEST(ParsingTest, test_Comment) {
  unsigned char file_content[] = {0x00, 0x13, 0x43, 0x72, 0x65, 0x61, 0x74,
                                  0x65, 0x64, 0x20, 0x77, 0x69, 0x74, 0x68,
                                  0x20, 0x47, 0x49, 0x4D, 0x50};
  unsigned int index = 0;
  std::string comment = ParseComment(file_content, index),
              expected("Created with GIMP");
  ASSERT_STREQ(comment.c_str(), expected.c_str());
  ASSERT_EQ(index, 19);
}

TEST(ParsingTest, test_FrameHeader) {
  unsigned char file_content[] = {0x00, 0x11, 0x08, 0x00, 0x10, 0x00,
                                  0x10, 0x03, 0x01, 0x11, 0x00, 0x02,
                                  0x11, 0x01, 0x03, 0x11, 0x01};
  unsigned int index = 0;
  std::vector<unsigned char> expected_results = {1, 1, 0},
                             expected_results_2 = {1, 1, 1},
                             expected_results_3 = {1, 1, 1};
  FrameHeader header =
      ParseFrameHeader(file_content, index, FRAME_TYPE_BASELINE_DTC);

  ASSERT_EQ(header.encoding_process_type_, FRAME_TYPE_BASELINE_DTC);
  ASSERT_EQ(header.sample_precision_, 8);
  ASSERT_EQ(header.number_of_lines_, 16);
  ASSERT_EQ(header.number_of_samples_per_line_, 16);
  ASSERT_EQ(header.number_of_component_, 3);
  ASSERT_THAT(header.component_parameters_.at(1),
              ::testing::ElementsAreArray(expected_results));
  ASSERT_THAT(header.component_parameters_.at(2),
              ::testing::ElementsAreArray(expected_results_2));
  ASSERT_THAT(header.component_parameters_.at(3),
              ::testing::ElementsAreArray(expected_results_3));
  ASSERT_EQ(index, 17);
}

TEST(ParsingTest, test_ScanHeader) {
  unsigned char file_content[] = {0x00, 0x0C, 0x03, 0x01, 0x00, 0x02,
                                  0x11, 0x03, 0x11, 0x00, 0x3F, 0x00};

  unsigned int index = 0;

  ScanHeader header = ParseScanHeader(file_content, index);

  ASSERT_EQ(header.number_of_component_s_, 3);
  ASSERT_EQ(header.start_of_spectral_selection_, 0);
  ASSERT_EQ(header.end_of_spectral_selection_, 63);
  ASSERT_EQ(header.approximation_high_bit_, 0);
  ASSERT_EQ(header.approximation_low_bit_, 0);

  ASSERT_EQ(header.components_parameters_.at(1).first, 0);
  ASSERT_EQ(header.components_parameters_.at(1).second, 0);
  ASSERT_EQ(header.components_parameters_.at(2).first, 1);
  ASSERT_EQ(header.components_parameters_.at(2).second, 1);
  ASSERT_EQ(header.components_parameters_.at(3).first, 1);
  ASSERT_EQ(header.components_parameters_.at(3).second, 1);

  ASSERT_EQ(index, 12);
}

TEST(ParsingTest, test_QuantizationTable) {
  unsigned char file_content[] = {0x00, 0x43, 0x01, 0x03, 0x04, 0x04, 0x05,
                                  0x04, 0x05, 0x09, 0x05, 0x05, 0x09, 0x14,
                                  0x0D, 0x0B, 0x0D, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
                                  0x14, 0x14, 0x14, 0x14},
                tq;
  unsigned int index = 0;
  std::vector<unsigned short> expected_qks = {
      0x03, 0x04, 0x04, 0x05, 0x04, 0x05, 0x09, 0x05, 0x05, 0x09, 0x14,
      0x0D, 0x0B, 0x0D, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
      0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
      0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
      0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
      0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14};

  QuantizationTable table;
  std::tie(tq, table) = ParseQuantizationTable(file_content, index);

  ASSERT_EQ(tq, 1);
  ASSERT_EQ(table.pq_, 0);
  ASSERT_THAT(table.qks_, ::testing::ElementsAreArray(expected_qks));
  ASSERT_EQ(index, 67);
}

TEST(ParsingTest, test_HuffmanTable) {
  unsigned char file_content[] = {0x00, 0x14, 0x10, 0x01, 0x00, 0x00, 0x00,
                                  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                  0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                table_index;
  unsigned int index = 0;
  std::vector<unsigned char> expected_bits = {0x01, 0x00, 0x00, 0x00,
                                              0x00, 0x00, 0x00, 0x00,
                                              0x00, 0x00, 0x00, 0x00,
                                              0x00, 0x00, 0x00, 0x00},
                             expected_huffvals = {0x00};
  HuffmanTable table;

  std::vector<std::pair<unsigned char, HuffmanTable>> results =
      ParseHuffmanTableSpecification(file_content, index);

  std::tie(table_index, table) = results.at(0);

  ASSERT_EQ(index, 20);
  ASSERT_EQ(table_index, 0);
  ASSERT_EQ(table.table_class_, 1);
  ASSERT_THAT(table.bits, ::testing::ElementsAreArray(expected_bits));
  ASSERT_THAT(table.huffvals, ::testing::ElementsAreArray(expected_huffvals));
}