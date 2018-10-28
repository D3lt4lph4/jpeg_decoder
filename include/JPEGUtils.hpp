#ifndef __JPEG_IMAGE__
#define __JPEG_IMAGE__

/**
 * \file JPEGUtils.hpp
 * \author Benjamin Deguerre
 * \brief File containing all the classes usefull for the decoder.
 */
#include <vector>

/**
 * \class JPEGImage
 * \brief A class representing an image decoded by a decoder. Depending on the
 * level of decoding, the components of an image may vary in size.
 */
class JPEGImage {
 public:
  JPEGImage();
  JPEGImage(std::vector<std::pair<int, int>> sizes);
  ~JPEGImage();

  std::pair<int, int> GetComponentShape(int component);
  std::vector<int> GetRealShape();
  int GetNumberOfComponent();
  void SetRealShape(std::vector<int> shape);
  void RescaleToRealSize();

  int& at(int row, int col, int component);
  std::vector<int>& GetData(int component);

 private:
  std::vector<std::vector<int>> image_components_;
  std::vector<std::pair<int, int>> components_shape;
  std::vector<int> real_shape_;
};

#endif