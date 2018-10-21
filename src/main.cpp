#include <boost/filesystem.hpp>
#include <cerrno>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <cxxopts.hpp>

#include "JPEGDecoder.hpp"

// We only use opencv if we are in debug mode to visually check the results of
// the decoding
#ifdef DEBUG
#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#endif

void matwrite(const std::string& filename, int* mat, unsigned int image_size_x,
              unsigned int image_size_y, int channels) {
  std::ofstream fs(filename, std::fstream::binary);
  if (fs.fail()) {
    std::cerr << strerror(errno) << std::endl;
  }

  fs.write((char*)&image_size_y, sizeof(int));  // rows
  fs.write((char*)&image_size_x, sizeof(int));  // cols
  fs.write((char*)&channels, sizeof(int));      // channels
  fs.write((char*)mat, sizeof(int) * image_size_x * image_size_y);
  fs.close();
}

int main(int argc, char* argv[]) {
#ifdef DEBUG
  bool show = false;
#endif
  unsigned int image_size_x, image_size_y;
  int channels;

  // Creating the parser.
  cxxopts::Options options("jpeg_decoder",
                           "Decoder for jpeg files, can decode at different "
                           "\"level\" and export the decoded image.");
  options.add_options()(
      "d,directory", "Directory containing the images to decode.",
      cxxopts::value<std::string>()->default_value("default"))(
      "l,level",
      "The level of decoding to use, one or below to output the huffman "
      "decoded values, 2 for the dequantized values, 3 for IDCT values, and 4 "
      "or more for the RGB values.",
      cxxopts::value<int>()->default_value("2"))(
      "f,file",
      "The file to parse, if a directory is specified, will be ignored.",
      cxxopts::value<std::string>()->default_value("default"))("help",
                                                               "Print help")
#ifdef DEBUG
      ("s,show",
       "If the image(s) should be displayed. For now, not handled by the "
       "directory parsing.",
       cxxopts::value<bool>(show))
#endif
          ("o,output",
           "The output directory, the generated file name will use the name of "
           "the "
           "image + .dat. If not specified, the program will try to output in "
           "the "
           "image directory.",
           cxxopts::value<std::string>());

  auto result = options.parse(argc, argv);

  // If the help is required, display it.
  if (result.count("help")) {
    std::cout << options.help({"", "Group"}) << std::endl;
    exit(0);
  }
  JPEGDecoder decoder;
  int* image;
  image = (int*)decoder.DecodeFile("data/chat_2.jpg", 4);
  if (!result.count("level")) {
    std::cout
        << "The level should be specified for the decoder to know when to stop."
        << std::endl;
    exit(0);
  }

  // If the directory is specified, we process the files in it.
  if (result.count("directory")) {
    std::string directory = result["directory"].as<std::string>();
    boost::filesystem::directory_iterator end_iterator;
    if (!boost::filesystem::exists(directory)) {
      std::cout << "The directory doesn't exists, exiting..." << std::endl;
      exit(0);
    }

    for (boost::filesystem::directory_iterator iterator(directory);
         iterator != end_iterator; ++iterator) {
      if (!boost::filesystem::is_directory(iterator->status())) {
        if (boost::filesystem::extension(iterator->path()) == ".jpg") {
          JPEGDecoder decoder;
          int* image;
          unsigned int image_size_x = 0, image_size_y = 0;

          std::cout << "Processing the image : " << iterator->path().string()
                    << std::endl;
          image = (int*)decoder.DecodeFile(iterator->path().string(),
                                           result["level"].as<int>());

          image_size_x = decoder.getImageSizeX();
          image_size_y = decoder.getImageSizeY();
          channels = decoder.getChannels();

          // Writing the decoded image as .dat file.
          if (result.count("output")) {
            switch (result["level"].as<int>()) {
              case 2:
                matwrite(result["output"].as<std::string>() +
                             iterator->path().stem().string() + ".qhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
              case 3:
                matwrite(result["output"].as<std::string>() +
                             iterator->path().stem().string() + ".iqhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
              default:
                matwrite(result["output"].as<std::string>() +
                             iterator->path().stem().string() + ".riqhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
            }
          } else {
            switch (result["level"].as<int>()) {
              case 2:
                matwrite(iterator->path().parent_path().string() +
                             iterator->path().stem().string() + ".qhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
              case 3:
                matwrite(iterator->path().parent_path().string() +
                             iterator->path().stem().string() + ".iqhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
              default:
                matwrite(iterator->path().parent_path().string() +
                             iterator->path().stem().string() + ".riqhjpg",
                         image, image_size_x, image_size_y, channels);
                break;
            }
          }
        }
      }
    }
    return 0;
  } else {
    std::cout << "No directory specified, defaulting to files." << std::endl;
  }

  // Process the file if specified.
  if (result.count("file")) {
    JPEGDecoder decoder;
    int* image;
    image = (int*)decoder.DecodeFile(result["file"].as<std::string>(),
                                     result["level"].as<int>());
    image_size_x = decoder.getImageSizeX();
    image_size_y = decoder.getImageSizeY();
    channels = decoder.getChannels();
    // Writing the decoded image as .dat file.
    boost::filesystem::path p(result["file"].as<std::string>());

    if (result.count("output")) {
      switch (result["level"].as<int>()) {
        case 2:
          matwrite(
              result["output"].as<std::string>() + p.stem().string() + ".qhjpg",
              image, image_size_x, image_size_y, channels);
          break;
        case 3:
          matwrite(result["output"].as<std::string>() + p.stem().string() +
                       ".iqhjpg",
                   image, image_size_x, image_size_y, channels);
          break;
        default:
          matwrite(result["output"].as<std::string>() + p.stem().string() +
                       ".riqhjpg",
                   image, image_size_x, image_size_y, channels);
          break;
      }
    } else {
      switch (result["level"].as<int>()) {
        case 2:
          matwrite(p.parent_path().string() + p.stem().string() + ".qhjpg",
                   image, image_size_x, image_size_y, channels);
          break;
        case 3:
          matwrite(p.parent_path().string() + p.stem().string() + ".iqhjpg",
                   image, image_size_x, image_size_y, channels);
          break;
        default:
          matwrite(p.parent_path().string() + p.stem().string() + ".riqhjpg",
                   image, image_size_x, image_size_y, channels);
          break;
      }
    }
#ifdef DEBUG
    if (show) {
      cv::Mat image_to_display =
          cv::Mat(decoder.getBlockPerColumn() * 8,
                  (decoder.getBlockPerLine()) * 8, CV_32SC3);
      int row_size = decoder.getBlockPerLine();
      std::cout << decoder.getBlockPerColumn() * 8 << std::endl;
      std::cout << (decoder.getBlockPerLine())* 8 << std::endl;
      image_size_y = 566;
      // Putting the data inside the opencv matrix.
      for (size_t row = 0; row < decoder.getBlockPerColumn(); row++) {
        for (size_t col = 0; col < decoder.getBlockPerLine(); col++) {
          for (size_t row_cell = 0; row_cell < 8; row_cell++) {
            for (size_t col_cell = 0; col_cell < 8; col_cell++) {
              image_to_display.at<cv::Vec3i>(row * 8 + row_cell,
                                             col * 8 + col_cell)[0] =
                  image[row * 64 * row_size * 3 + col * 64 * 3 + row_cell * 8 +
                        col_cell];
              image_to_display.at<cv::Vec3i>(row * 8 + row_cell,
                                             col * 8 + col_cell)[1] =
                  image[row * 64 * row_size * 3 + col * 64 * 3 + row_cell * 8 +
                        col_cell + 64];
              image_to_display.at<cv::Vec3i>(row * 8 + row_cell,
                                             col * 8 + col_cell)[2] =
                  image[row * 64 * row_size * 3 + col * 64 * 3 + row_cell * 8 +
                        col_cell + 128];
            }
          }
        }
      }

      image_to_display.convertTo(image_to_display, CV_8UC3);
      cv::imshow("Decoded image.", image_to_display);
      cv::waitKey(0);
    }
#endif
    return 0;
  } else {
    std::cout << "No file specified, exiting." << std::endl;
  }

  return 0;
}
