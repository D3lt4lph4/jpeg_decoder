#include "JPEGUtils.hpp"

JPEGImage::JPEGImage() {}

JPEGImage::JPEGImage(std::vector<std::pair<int, int>> sizes) {
  this->image_components_.resize(sizes.size());
  this->components_shape = sizes;

  for (int component = 0; component < sizes.size(); component++) {
    this->image_components_.at(component) =
        new int[sizes.at(component).first * sizes.at(component).second];
  }
}

JPEGImage::~JPEGImage() {
  for (int* object : this->image_components_) {
    delete[] object;
  }
  this->image_components_.clear();
}

std::pair<int, int> JPEGImage::GetComponentShape(int component) {
  return this->components_shape.at(component);
}

int& JPEGImage::at(int row, int col, int component) {
  int line_length = this->components_shape.at(component).first;
  return this->image_components_.at(component)[row * line_length + col];
}
