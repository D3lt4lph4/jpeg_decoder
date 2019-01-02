/**
 * \file JPEGUtils.cpp
 * \author Deguerre Benjamin
 * \brief Contains all the object to be used by the JPEGDecoder.
 */

#include "JPEGUtils.hpp"

/**
 * \class JPEGImage::JPEGImage
 * \brief A class representing an image decoded by a decoder. Depending on the
 * level of decoding, the components of an image may vary in size.
 */

/**
 * \fn JPEGImage()
 * \brief Dummy constructor for a JPEGImage.
 */
JPEGImage::JPEGImage() {}

/**
 * \fn JPEGImage::JPEGImage(std::vector<std::pair<int, int>> sizes)
 *
 * \brief Constructor of a JPEGImage, resize the vector with the correct sizes
 * for the data.
 *
 * \param[in] std::vector<std::pair<int, int>> sizes A vector containing the
 * sizes of each components, should be as (nb_components, row, col).
 */
JPEGImage::JPEGImage(std::vector<std::pair<int, int>> sizes) {
  this->image_components_.resize(sizes.size());
  this->components_shape = sizes;

  for (int component = 0; component < sizes.size(); component++) {
#ifdef DEBUG
    this->image_components_.at(component).resize(sizes.at(component).first *
                                                 sizes.at(component).second);
#else
    this->image_components_.at(component).resize(sizes.at(component).first *
                                                 sizes.at(component).second);
#endif
  }
}

/**
 * \fn JPEGImage::~JPEGImage()
 * \brief The destructor for the Image.
 */
JPEGImage::~JPEGImage() {}

/**
 * \fn std::pair<int, int> JPEGImage::GetComponentShape(int component)
 *
 * \brief Return the shape of the component specified. This is the shape of the
 component as seen by the decoder, i.e with the block padding and the
 subsampling.
 *
 * \param component The position of the component.
 *
 * \return The shape of the data at component, in the order (row, col).
 */
std::pair<int, int> JPEGImage::GetComponentShape(int component) {
#ifdef DEBUG
  return this->components_shape.at(component);
#else
  return this->components_shape.at(component);
#endif
}

/**
 * \fn JPEGImage::GetRealShape()
 * \brief Return the real shape of the image as it was before going through the
 * steps of subsampling.
 *
 * \return A vector containing the shape of the image. In the order, (row, col,
 * component).
 */
std::vector<int> JPEGImage::GetRealShape() { return this->real_shape_; }

/**
 * \fn int JPEGImage::GetNumberOfComponent()
 *
 * \return The number of components in the image.
 */
int JPEGImage::GetNumberOfComponent() { return this->image_components_.size(); }

/**
 * \fn void JPEGImage::SetRealShape(std::vector<int> shape)
 * \brief Setter for the real shape of the image.
 *
 * \param[in] shape A vector of int containing the real shape of the image (row,
 * col, nb_components)
 */
void JPEGImage::SetRealShape(std::vector<int> shape) {
  this->real_shape_ = shape;
}

/**
 * \fn int& JPEGImage::at(int row, int col, int component)
 * \brief Return the value at the place indicated. The value is return by
 * reference, i.e alterable.
 *
 * \param[in] row The position of the row
 * \param[in] col The position of the column
 * \param[in] component The position of the component
 *
 * \return The value at (row, col, component)
 */
int& JPEGImage::at(int row, int col, int component) {
#ifdef DEBUG
  int line_length = this->components_shape.at(component).second;
  return this->image_components_.at(component).at(row * line_length + col);
#else
  int line_length = this->components_shape.at(component).second;
  return this->image_components_.at(component).at(row * line_length + col);
#endif
}

/**
 * \fn std::vector<int>& JPEGImage::GetData(int component)
 * \brief Return the vector<int> representing the data at component. The data is
 stored line by line in a vector.
 *
 * \return The vector containing the data, if the returned value is modified,
 the vector in the image will be too.
 */
std::vector<int>& JPEGImage::GetData(int component) {
#ifdef DEBUG
  return this->image_components_.at(component);
#else
  return this->image_components_.at(component);
#endif
}

/**
 * \fn void JPEGImage::RescaleToRealSize()
 * \brief Rescale the image in the Object to the real size. This function will
 * resize each component if required. For now only the upsizing rezise to the
 * first components.
 */
void JPEGImage::RescaleToRealSize() {
  int col_factor, row_factor;

  int rows = this->components_shape.at(0).first,
      cols = this->components_shape.at(0).second;
  for (int i = 1; i < this->components_shape.size(); i++) {
    std::vector<int> new_data(rows * cols);
    col_factor = cols / this->components_shape.at(i).second;
    row_factor = rows / this->components_shape.at(i).first;
    for (int row = 0; row < this->components_shape.at(i).first; row++) {
      for (int col = 0; col < this->components_shape.at(i).second; col++) {
        for (int row_f = 0; row_f < row_factor; row_f++) {
          for (int col_f = 0; col_f < col_factor; col_f++) {
#ifdef DEBUG
            new_data.at((row * row_factor + row_f) * cols + col * col_factor +
                        col_f) =
                this->image_components_.at(i).at(
                    row * this->components_shape.at(i).second + col);
#else
            new_data[(row * row_factor + row_f) * cols + col * col_factor +
                     col_f] =
                this->image_components_[i]
                                       [row * this->components_shape[i].second +
                                        col];
#endif
          }
        }
      }
    }
    this->image_components_.at(i) = new_data;
#ifdef DEBUG
    this->components_shape.at(i).first = rows;
    this->components_shape.at(i).second = cols;
#else
    this->components_shape[i].first = rows;
    this->components_shape[i].second = cols;
#endif
  }
}