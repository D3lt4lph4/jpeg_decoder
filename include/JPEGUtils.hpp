#ifndef __JPEG_IMAGE__
#define __JPEG_IMAGE__

#include <vector>
#include "JPEGType.hpp"

class JPEGImage {
 public:
  JPEGImage();
  JPEGImage(std::vector<std::pair<int, int>> sizes);
  ~JPEGImage();

  std::pair<int, int> GetComponentShape(int component);
  std::vector<int> GetRealShape();
  int GetNumberOfComponent();
  ImageType GetType();
  void SetType(ImageType type);
  void SetRealShape(std::vector<int> shape);
  void RescaleToRealSize();

  int& at(int row, int col, int component);
  std::vector<int>& GetData(int component);

 private:
  std::vector<std::vector<int>> image_components_;
  std::vector<std::pair<int, int>> components_shape;
  std::vector<int> real_shape_;
  ImageType type_;
};

#endif