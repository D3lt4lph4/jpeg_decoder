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

TEST(TestHuffmanFunctions, testDecodeACCoefficients) {}

TEST(TestHuffmanFunctions, testDecodeZZ) {}
