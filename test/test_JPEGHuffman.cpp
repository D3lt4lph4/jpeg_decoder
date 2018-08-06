#include "JPEGDecoder.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

class AdditionTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST(TestHuffmanFunctions, testGenerateSizeTable) {
  std::vector<unsigned char> bits = {1, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0},
                             bits_2 = {0, 1, 5, 1, 1, 1, 1, 1,
                                       1, 0, 0, 0, 0, 0, 0, 0},
                             size_table = {1, 0},
                             size_table_2 = {2, 3, 3, 3, 3, 3, 4,
                                             5, 6, 7, 8, 9, 0},
                             size_table_res, size_table_res_2;

  unsigned char last_k, last_k_2;

  std::tie(last_k, size_table_res) = GenerateSizeTable(bits);

  ASSERT_EQ(last_k, 1);
  ASSERT_THAT(size_table_res, ::testing::ElementsAreArray(size_table));

  std::tie(last_k_2, size_table_res_2) = GenerateSizeTable(bits_2);

  ASSERT_EQ(last_k_2, 12);
  ASSERT_THAT(size_table_res_2, ::testing::ElementsAreArray(size_table_2));
}

TEST(TestHuffmanFunctions, testGenerateCodeTable) {
  std::vector<unsigned char> huffsize = {1, 0},
                             huffsize_2 = {2, 3, 3, 3, 3, 3, 4,
                                           5, 6, 7, 8, 9, 0};

  std::vector<unsigned short> code_table = {0},
                              code_table_2 = {0,  2,  3,  4,   5,   6,
                                              14, 30, 62, 126, 254, 510},
                              code_table_res, code_table_res_2;

  code_table_res = GenerateCodeTable(huffsize);
  ASSERT_THAT(code_table_res, ::testing::ElementsAreArray(code_table));

  code_table_res_2 = GenerateCodeTable(huffsize_2);
  ASSERT_THAT(code_table_res_2, ::testing::ElementsAreArray(code_table_2));
}

TEST(TestHuffmanFunctions, testDecoderTable) {
  std::vector<unsigned short> code_table = {0},
                              code_table_2 = {0,  2,  3,  4,   5,   6,
                                              14, 30, 62, 126, 254, 510};

  std::vector<int> min_code = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0},
                   min_code_2 = {0,   0, 0, 2, 14, 30, 62, 126, 254,
                                 510, 0, 0, 0, 0,  0,  0,  0},
                   max_code = {0,  0,  -1, -1, -1, -1, -1, -1, -1,
                               -1, -1, -1, -1, -1, -1, -1, -1},
                   max_code_2 = {0,   -1, 0,  6,  14, 30, 62, 126, 254,
                                 510, -1, -1, -1, -1, -1, -1, -1},
                   min_code_res, min_code_2_res, max_code_res, max_code_2_res;

  std::vector<unsigned char> bits = {1, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0},
                             bits_2 = {0, 1, 5, 1, 1, 1, 1, 1,
                                       1, 0, 0, 0, 0, 0, 0, 0},
                             valptr = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0},
                             valptr_2 = {0,  0, 0, 1, 6, 7, 8, 9, 10,
                                         11, 0, 0, 0, 0, 0, 0, 0},
                             valptr_res, valptr_2_res;

  std::tie(max_code_res, min_code_res, valptr_res) =
      DecoderTables(bits, code_table);
  std::tie(max_code_2_res, min_code_2_res, valptr_2_res) =
      DecoderTables(bits_2, code_table_2);

  std::cout << "testing for the first set of elements" << std::endl;
  ASSERT_THAT(min_code_res, ::testing::ElementsAreArray(min_code));
  ASSERT_THAT(max_code_res, ::testing::ElementsAreArray(max_code));
  ASSERT_THAT(valptr_res, ::testing::ElementsAreArray(valptr));

  ASSERT_THAT(min_code_2_res, ::testing::ElementsAreArray(min_code_2));
  ASSERT_THAT(max_code_2_res, ::testing::ElementsAreArray(max_code_2));
  ASSERT_THAT(valptr_2_res, ::testing::ElementsAreArray(valptr_2));
}

TEST(TestHuffmanFunctions, testDecode) {
  HuffmanTable test_table;
  unsigned char res, file_content[4] = {0x6E, 0x43, 0xF0, 0x65}, bit_index = 8;
  unsigned int index = 0;
  std::vector<unsigned char> valptr = {0,  0, 0, 1, 6, 7, 8, 9, 10,
                                       11, 0, 0, 0, 0, 0, 0, 0},
                             bits = {0, 1, 5, 1, 1, 1, 1, 1,
                                     1, 0, 0, 0, 0, 0, 0, 0},
                             size_table = {2, 3, 3, 3, 3, 3, 4,
                                           5, 6, 7, 8, 9, 0},
                             huffvals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  std::vector<unsigned short> code_table = {0,  2,  3,  4,   5,   6,
                                            14, 30, 62, 126, 254, 510};
  std::vector<int> min_code = {0,   0, 0, 2, 14, 30, 62, 126, 254,
                               510, 0, 0, 0, 0,  0,  0,  0},
                   max_code = {0,   -1, 0,  6,  14, 30, 62, 126, 254,
                               510, -1, -1, -1, -1, -1, -1, -1};
  test_table.bits = bits;
  test_table.huffsize = size_table;
  test_table.huffvals = huffvals;
  test_table.huffcode = code_table;
  test_table.min_code = min_code;
  test_table.max_code = max_code;
  test_table.val_pointer = valptr;

  res = Decode(file_content, &index, &bit_index, test_table);
  ASSERT_EQ(res, 2);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 0);

  bit_index = 4;
  res = Decode(file_content, &index, &bit_index, test_table);
  ASSERT_EQ(res, 6);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 1);

  res = Decode(file_content, &index, &bit_index, test_table);
  ASSERT_EQ(res, 1);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 1);

  res = Decode(file_content, &index, &bit_index, test_table);
  ASSERT_EQ(res, 0);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 1);

  bit_index = 2;
  res = Decode(file_content, &index, &bit_index, test_table);
  ASSERT_EQ(res, 9);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 2);
}

TEST(TestHuffmanFunctions, testReceive) {
  unsigned char file_content[4] = {0x6E, 0x43, 0xF0, 0x65}, bit_index = 7;
  int res;
  unsigned int index = 0;

  res = Receive(9, file_content, &index, &bit_index);

  ASSERT_EQ(res, 441);
  ASSERT_EQ(bit_index, 6);
  ASSERT_EQ(index, 1);

  res = Receive(3, file_content, &index, &bit_index);

  ASSERT_EQ(res, 0);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 1);

  res = Receive(3, file_content, &index, &bit_index);

  ASSERT_EQ(res, 3);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 2);
}

TEST(TestHuffmanFunctions, testExtended) {
  ASSERT_EQ(Extended(441, 9), 441);
  ASSERT_EQ(Extended(185, 9), -326);
  ASSERT_EQ(Extended(185, 8), 185);
  ASSERT_EQ(Extended(57, 7), -70);
  ASSERT_EQ(Extended(57, 6), 57);
  ASSERT_EQ(Extended(1, 1), 1);
  ASSERT_EQ(Extended(0, 1), -1);
}

TEST(TestHuffmanFunctions, testDecodeACCoefficients) {
  unsigned char file_content[6] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                bit_index = 8;
  unsigned int index = 0;

  std::vector<unsigned char>
      valptr = {0, 0, 0, 1, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0, 0, 0},
      bits = {0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
      bits_2 = {0x00, 0x02, 0x01, 0x02, 0x04, 0x04, 0x03, 0x04,
                0x07, 0x05, 0x04, 0x04, 0x00, 0x01, 0x02, 0x77},
      size_table = {2, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 0},
      huffvals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
      huffvals_2 = {
          0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12,
          0x41, 0x51, 0x07, 0x61, 0x71, 0x13, 0x22, 0x32, 0x81, 0x08, 0x14,
          0x42, 0x91, 0xA1, 0xB1, 0xC1, 0x09, 0x23, 0x33, 0x52, 0xF0, 0x15,
          0x62, 0x72, 0xD1, 0x0A, 0x16, 0x24, 0x34, 0xE1, 0x25, 0xF1, 0x17,
          0x18, 0x19, 0x1A, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x35, 0x36, 0x37,
          0x38, 0x39, 0x3A, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4A,
          0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x63, 0x64, 0x65,
          0x66, 0x67, 0x68, 0x69, 0x6A, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
          0x79, 0x7A, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8A,
          0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0xA2, 0xA3,
          0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xB2, 0xB3, 0xB4, 0xB5,
          0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xC2, 0xC3, 0xC4, 0xC5, 0xC6, 0xC7,
          0xC8, 0xC9, 0xCA, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9,
          0xDA, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xF2,
          0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8, 0xF9, 0xFA};
  std::vector<unsigned short> code_table = {0,  2,  3,  4,   5,   6,
                                            14, 30, 62, 126, 254, 510};
  std::vector<int> min_code = {0,   0, 0, 2, 14, 30, 62, 126, 254,
                               510, 0, 0, 0, 0,  0,  0,  0},
                   max_code = {0,   -1, 0,  6,  14, 30, 62, 126, 254,
                               510, -1, -1, -1, -1, -1, -1, -1},
                   results,
                   expected_results = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  HuffmanTable first_table, second_table;

  first_table.bits = bits;
  first_table.huffsize = size_table;
  first_table.huffvals = huffvals;
  first_table.huffcode = code_table;
  first_table.min_code = min_code;
  first_table.max_code = max_code;
  first_table.val_pointer = valptr;

  second_table.bits = bits_2;
  std::cout << "azer" << std::endl;
  std::tie(std::ignore, second_table.huffsize) = GenerateSizeTable(bits_2);
  std::cout << "azer" << std::endl;
  second_table.huffvals = huffvals_2;
  std::cout << "azer" << std::endl;
  second_table.huffcode = GenerateCodeTable(second_table.huffsize);
  std::cout << "qsdf" << std::endl;
  std::tie(second_table.max_code, second_table.min_code,
           second_table.val_pointer) =
      DecoderTables(bits_2, second_table.huffcode);
  std::cout << "qsdf" << std::endl;
  // First let's test the decoding of a full 0's block.
  results = DecodeACCoefficients(file_content, &index, &bit_index, first_table);
  std::cout << "ret" << std::endl;
  ASSERT_THAT(results, ::testing::ElementsAreArray(expected_results));
  ASSERT_EQ(index, 0);
  ASSERT_EQ(bit_index, 6);
  std::cout << "rty" << std::endl;
  index = 0;
  bit_index = 8;
  results =
      DecodeACCoefficients(file_content, &index, &bit_index, second_table);
  std::cout << "azer" << std::endl;
  ASSERT_THAT(results, ::testing::ElementsAreArray(expected_results));

  // Then we test with a block of 4 black pixels followed by 4 white, 4 black,
  // ...
}

TEST(TestHuffmanFunctions, testDecodeZZ) {
  unsigned char file_content[4] = {0xDC, 0xDC, 0xAA, 0xAA}, bit_index = 8;
  int res;
  unsigned int index = 0;

  ASSERT_EQ(DecodeZZ(file_content, &index, &bit_index, 9), 441);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 1);

  ASSERT_EQ(DecodeZZ(file_content, &index, &bit_index, 8), 185);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 2);
}
