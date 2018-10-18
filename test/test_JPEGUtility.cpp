#include "JPEGUtility.hpp"
#include "gmock/gmock.h"

TEST(UtilityTest, gettingNextBit) {
  unsigned char bit, file_content, bit_index;
  unsigned int index;

  file_content = 0x93;
  bit_index = 8;
  index = 0;

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 6);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 4);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 2);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 1);
  ASSERT_EQ(index, 0);

  bit = NextBit(&file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 1);
}

TEST(UtilityTest, incrementOnFF00) {
  unsigned char bit, file_content[3] = {0xFF, 0x00, 0x23}, bit_index = 8;
  unsigned int index = 0;

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 6);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 5);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 4);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 3);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 2);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 1);
  ASSERT_EQ(index, 0);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 1);
  ASSERT_EQ(bit_index, 8);
  ASSERT_EQ(index, 2);

  bit = NextBit(file_content, &index, &bit_index);
  ASSERT_EQ(bit, 0);
  ASSERT_EQ(bit_index, 7);
  ASSERT_EQ(index, 2);
}

TEST(UtilityTest, testThrow) {
  unsigned char bit, file_content[3] = {0xFF, 0x34, 0x65}, bit_index = 0;
  unsigned int index = 0;

  EXPECT_THROW(NextBit(file_content, &index, &bit_index), std::out_of_range);

  bit_index = 9;
  EXPECT_THROW(NextBit(file_content, &index, &bit_index), std::out_of_range);

  bit_index = 8;
  EXPECT_THROW(NextBit(file_content, &index, &bit_index), std::runtime_error);
}

TEST(UtilityTest, testIDCT) {
  int value_input[] = {0,  0,  0,  92, 92, 92, 0,  0,  0,  72, 72, 72, 0, 0, 0,
                       41, 41, 41, 0,  0,  0,  77, 77, 77, 0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},

      value_input_2[] = {
          -4,  -4,  -4, 184, 184, 184, 0,   0,   0, 215, 215, 215, 0, 0, 0, 326,
          326, 326, 0,  0,   0,   921, 921, 921, 0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0};

  int value_result[] = {161, 161, 161, 124, 124, 124, 137, 137, 137, 116, 116,
                        116, 139, 139, 139, 118, 118, 118, 131, 131, 131, 94,
                        94,  94,  161, 161, 161, 124, 124, 124, 137, 137, 137,
                        116, 116, 116, 139, 139, 139, 118, 118, 118, 131, 131,
                        131, 94,  94,  94,  161, 161, 161, 124, 124, 124, 137,
                        137, 137, 116, 116, 116, 139, 139, 139, 118, 118, 118,
                        131, 131, 131, 94,  94,  94,  161, 161, 161, 124, 124,
                        124, 137, 137, 137, 116, 116, 116, 139, 139, 139, 118,
                        118, 118, 131, 131, 131, 94,  94,  94,  161, 161, 161,
                        124, 124, 124, 137, 137, 137, 116, 116, 116, 139, 139,
                        139, 118, 118, 118, 131, 131, 131, 94,  94,  94,  161,
                        161, 161, 124, 124, 124, 137, 137, 137, 116, 116, 116,
                        139, 139, 139, 118, 118, 118, 131, 131, 131, 94,  94,
                        94,  161, 161, 161, 124, 124, 124, 137, 137, 137, 116,
                        116, 116, 139, 139, 139, 118, 118, 118, 131, 131, 131,
                        94,  94,  94,  161, 161, 161, 124, 124, 124, 137, 137,
                        137, 116, 116, 116, 139, 139, 139, 118, 118, 118, 131,
                        131, 131, 94,  94,  94},
      value_result_2[] = {
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0};

  IDCT(value_input);
  IDCT(value_input+64);
  IDCT(value_input+128);
  IDCT(value_input_2);
  IDCT(value_input_2+64);
  IDCT(value_input_2+128);

  ASSERT_THAT(value_input, ::testing::ElementsAreArray(value_result));
  ASSERT_THAT(value_input_2, ::testing::ElementsAreArray(value_result_2));
}

TEST(UtilityTest, testFastIDCT) {
  int value_input[] = {0,  0,  0,  92, 92, 92, 0,  0,  0,  72, 72, 72, 0, 0, 0,
                       41, 41, 41, 0,  0,  0,  77, 77, 77, 0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},

      value_input_2[] = {
          -4,  -4,  -4, 184, 184, 184, 0,   0,   0, 215, 215, 215, 0, 0, 0, 326,
          326, 326, 0,  0,   0,   921, 921, 921, 0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0,
          0,   0,   0,  0,   0,   0,   0,   0,   0, 0,   0,   0,   0, 0, 0, 0};

  int value_result[] = {161, 161, 161, 124, 124, 124, 137, 137, 137, 116, 116,
                        116, 139, 139, 139, 118, 118, 118, 131, 131, 131, 94,
                        94,  94,  161, 161, 161, 124, 124, 124, 137, 137, 137,
                        116, 116, 116, 139, 139, 139, 118, 118, 118, 131, 131,
                        131, 94,  94,  94,  161, 161, 161, 124, 124, 124, 137,
                        137, 137, 116, 116, 116, 139, 139, 139, 118, 118, 118,
                        131, 131, 131, 94,  94,  94,  161, 161, 161, 124, 124,
                        124, 137, 137, 137, 116, 116, 116, 139, 139, 139, 118,
                        118, 118, 131, 131, 131, 94,  94,  94,  161, 161, 161,
                        124, 124, 124, 137, 137, 137, 116, 116, 116, 139, 139,
                        139, 118, 118, 118, 131, 131, 131, 94,  94,  94,  161,
                        161, 161, 124, 124, 124, 137, 137, 137, 116, 116, 116,
                        139, 139, 139, 118, 118, 118, 131, 131, 131, 94,  94,
                        94,  161, 161, 161, 124, 124, 124, 137, 137, 137, 116,
                        116, 116, 139, 139, 139, 118, 118, 118, 131, 131, 131,
                        94,  94,  94,  161, 161, 161, 124, 124, 124, 137, 137,
                        137, 116, 116, 116, 139, 139, 139, 118, 118, 118, 131,
                        131, 131, 94,  94,  94},
      value_result_2[] = {
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,
          0,   0,   254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,
          254, 254, 254, 0,   0,   0,   254, 254, 254, 0,   0,   0,   254, 254,
          254, 0,   0,   0,   254, 254, 254, 0,   0,   0};
  
  
  FastIDCT(value_input);
  FastIDCT(value_input+64);
  FastIDCT(value_input+128);
  FastIDCT(value_input_2);
  FastIDCT(value_input_2+64);
  FastIDCT(value_input_2+128);

  ASSERT_THAT(value_input, ::testing::ElementsAreArray(value_result));
  ASSERT_THAT(value_input_2, ::testing::ElementsAreArray(value_result_2));
}

TEST(UtilityTest, testFastIDCT1) {

  int value_input[] = {-4, 184, 0, 215, 0, 326, 0, 921};

  int value_result[] = {358.60, -361.60, 358.96, -359.31, 356.48, -361.78,  358.77, -361.43};

  int results[8];

  FastIDCT1(value_input, results, 9, 1 << 8);
  
  ASSERT_THAT(results, ::testing::ElementsAreArray(value_result));
}