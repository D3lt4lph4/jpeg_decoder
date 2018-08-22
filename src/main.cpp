#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGDecoder.hpp"

void matwrite(const std::string& filename, const cv::Mat& mat) {
  std::ofstream fs(filename, std::fstream::binary);

  // Header
  int type = mat.type();
  int channels = mat.channels();
  fs.write((char*)&mat.rows, sizeof(int));  // rows
  fs.write((char*)&mat.cols, sizeof(int));  // cols
  fs.write((char*)&type, sizeof(int));      // type
  fs.write((char*)&channels, sizeof(int));  // channels

  // Data
  if (mat.isContinuous()) {
    fs.write(mat.ptr<char>(0), (mat.dataend - mat.datastart));
  } else {
    int row_size = CV_ELEM_SIZE(type) * mat.cols;
    for (int r = 0; r < mat.rows; ++r) {
      fs.write(mat.ptr<char>(r), row_size);
    }
  }
}

int main(int argc, char const* argv[]) {
  if (argc == 2) {
    JPEGDecoder decoder;
    cv::Mat image;
    std::string file_name = argv[1];

    image = decoder.DecodeFile(file_name, 2);
    matwrite("output.dat", image);
    image.convertTo(image, CV_8UC3);

    cv::imshow("Decoded image.", image);
    cv::waitKey(0);

  } else {
    std::cout << "Please enter one image to decode." << std::endl;
  }

  return 0;
}
