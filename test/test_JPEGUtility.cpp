#include "JPEGUtility.hpp"
#include "gtest/gtest.h"

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