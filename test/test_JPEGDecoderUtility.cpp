#include "JPEGDecoder.hpp"
#include "gtest/gtest.h"

class AdditionTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(AdditionTest, twoValues) {
  const int x = 4;
  const int y = 5;
}