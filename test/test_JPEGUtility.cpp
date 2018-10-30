#include "JPEGUtility.hpp"
#include "gmock/gmock.h"

TEST(UtilityTest, gettingNextBit) {
  unsigned char bit, file_content, bit_index;
  unsigned int index;

  file_content = 0x93;
  bit_index = 8;
  index = 0;

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 6);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 4);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 2);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 1);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 1);
}

TEST(UtilityTest, incrementOnFF00) {
  unsigned char bit, file_content[3] = {0xFF, 0x00, 0x23}, bit_index = 8;
  unsigned int index = 0;

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 6);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 4);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 2);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 1);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 2);

  bit = NextBit(file_content, index, bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 2);
}

TEST(UtilityTest, testThrow) {
  unsigned char bit, file_content[3] = {0xFF, 0x34, 0x65}, bit_index = 0;
  unsigned int index = 0;

  EXPECT_THROW(NextBit(file_content, index, bit_index), std::out_of_range);

  bit_index = 9;
  EXPECT_THROW(NextBit(file_content, index, bit_index), std::out_of_range);

  bit_index = 8;
  EXPECT_THROW(NextBit(file_content, index, bit_index), std::runtime_error);
}

TEST(UtilityTest, testIDCT) {
  int value_input[] = {0, 92, 0, 72, 0, 41, 0, 77, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 92, 0, 72, 0, 41, 0, 77, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 92, 0, 72, 0, 41, 0, 77, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0},

      value_input_2[] = {
          -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0};

  int value_result[] = {161, 124, 137, 116, 139, 118, 131, 94,  161, 124, 137,
                        116, 139, 118, 131, 94,  161, 124, 137, 116, 139, 118,
                        131, 94,  161, 124, 137, 116, 139, 118, 131, 94,  161,
                        124, 137, 116, 139, 118, 131, 94,  161, 124, 137, 116,
                        139, 118, 131, 94,  161, 124, 137, 116, 139, 118, 131,
                        94,  161, 124, 137, 116, 139, 118, 131, 94,  161, 124,
                        137, 116, 139, 118, 131, 94,  161, 124, 137, 116, 139,
                        118, 131, 94,  161, 124, 137, 116, 139, 118, 131, 94,
                        161, 124, 137, 116, 139, 118, 131, 94,  161, 124, 137,
                        116, 139, 118, 131, 94,  161, 124, 137, 116, 139, 118,
                        131, 94,  161, 124, 137, 116, 139, 118, 131, 94,  161,
                        124, 137, 116, 139, 118, 131, 94,  161, 124, 137, 116,
                        139, 118, 131, 94,  161, 124, 137, 116, 139, 118, 131,
                        94,  161, 124, 137, 116, 139, 118, 131, 94,  161, 124,
                        137, 116, 139, 118, 131, 94,  161, 124, 137, 116, 139,
                        118, 131, 94,  161, 124, 137, 116, 139, 118, 131, 94,
                        161, 124, 137, 116, 139, 118, 131, 94,  161, 124, 137,
                        116, 139, 118, 131, 94},
      value_result_2[] = {
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0,
          254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0, 254, 0};

  IDCT(value_input);
  IDCT(value_input + 64);
  IDCT(value_input + 128);
  IDCT(value_input_2);
  IDCT(value_input_2 + 64);
  IDCT(value_input_2 + 128);

  ASSERT_THAT(value_input, ::testing::ElementsAreArray(value_result));
  ASSERT_THAT(value_input_2, ::testing::ElementsAreArray(value_result_2));
}

TEST(UtilityTest, testFastIDCT) {
  std::vector<int> value_input = {0, 92, 0, 72, 0, 41, 0, 77, 0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 92, 0, 72, 0, 41, 0, 77,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 92, 0, 72,
                                  0, 41, 0, 77, 0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                                  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0},

                   value_input_2 = {
                       -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       -4, 184, 0, 215, 0, 326, 0, 921, 0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                       0,  0,   0, 0,   0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0};

  int value_result[] = {161, 124, 137, 117, 139, 119, 132, 95,  161, 124, 137,
                        117, 139, 119, 132, 95,  161, 124, 137, 117, 139, 119,
                        132, 95,  161, 124, 137, 117, 139, 119, 132, 95,  161,
                        124, 137, 117, 139, 119, 132, 95,  161, 124, 137, 117,
                        139, 119, 132, 95,  161, 124, 137, 117, 139, 119, 132,
                        95,  161, 124, 137, 117, 139, 119, 132, 95,  161, 124,
                        137, 117, 139, 119, 132, 95,  161, 124, 137, 117, 139,
                        119, 132, 95,  161, 124, 137, 117, 139, 119, 132, 95,
                        161, 124, 137, 117, 139, 119, 132, 95,  161, 124, 137,
                        117, 139, 119, 132, 95,  161, 124, 137, 117, 139, 119,
                        132, 95,  161, 124, 137, 117, 139, 119, 132, 95,  161,
                        124, 137, 117, 139, 119, 132, 95,  161, 124, 137, 117,
                        139, 119, 132, 95,  161, 124, 137, 117, 139, 119, 132,
                        95,  161, 124, 137, 117, 139, 119, 132, 95,  161, 124,
                        137, 117, 139, 119, 132, 95,  161, 124, 137, 117, 139,
                        119, 132, 95,  161, 124, 137, 117, 139, 119, 132, 95,
                        161, 124, 137, 117, 139, 119, 132, 95,  161, 124, 137,
                        117, 139, 119, 132, 95},
      value_result_2[] = {
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0,
          255, 0, 255, 1, 254, 0, 255, 0, 255, 0, 255, 1, 254, 0, 255, 0};

  FastIDCT2D(value_input, 0, 0, 8);
  FastIDCT2D(value_input, 8, 0, 8);
  FastIDCT2D(value_input, 16, 0, 8);
  FastIDCT2D(value_input_2, 0, 0, 8);
  FastIDCT2D(value_input_2, 8, 0, 8);
  FastIDCT2D(value_input_2, 16, 0, 8);

  ASSERT_THAT(value_input, ::testing::ElementsAreArray(value_result));
  ASSERT_THAT(value_input_2, ::testing::ElementsAreArray(value_result_2));
}

TEST(UtilityTest, testFastIDCT1) {
  std::vector<int> value_input = {-4, 184, 0, 215, 0, 326, 0, 921},
                   value_result = {1015, -1022, 1015, -1015,
                                   1007, -1023, 1014, -1023},
                   results(64);

  FastIDCT1D(value_input, results, 0, 0, 9, 1 << 8, 8);

  for (size_t i = 0; i < 8; i++) {
    results[i] = results[i * 8];
  }
  results.resize(8);

  ASSERT_THAT(results, ::testing::ElementsAreArray(value_result));
}