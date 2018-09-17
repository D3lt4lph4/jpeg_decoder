#include <boost/filesystem.hpp>
#include <cerrno>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <cxxopts.hpp>

#include "JPEGDecoder.hpp"

void matwrite(const std::string& filename, int* mat, unsigned int image_size_x,
              unsigned int image_size_y) {
  double min, max;
  std::ofstream fs(filename, std::fstream::binary);
  if (fs.fail()) {
    std::cerr << strerror(errno) << std::endl;
  }

  int channels = 3;
  fs.write((char*)image_size_y, sizeof(unsigned int));  // rows
  fs.write((char*)image_size_x, sizeof(int));           // cols
  fs.write((char*)&channels, sizeof(int));              // channels

  fs.write((char*)mat, image_size_x * image_size_y * sizeof(int));
}

int main(int argc, char* argv[]) {
  bool show = false;

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
      cxxopts::value<std::string>()->default_value("default"))(
      "help", "Print help")("s,show",
                            "If the image(s) should be displayed. For now, not "
                            "handled by the directory parsing.",
                            cxxopts::value<bool>(show))(
      "o,output",
      "The output directory, the generated file name will use the name of the "
      "image + .dat. If not specified, the program will try to output in the "
      "image directory.",
      cxxopts::value<std::string>());

  auto result = options.parse(argc, argv);

  // If the help is required, display it.
  if (result.count("help")) {
    std::cout << options.help({"", "Group"}) << std::endl;
    exit(0);
  }

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
                                           &image_size_x, &image_size_y,
                                           result["level"].as<int>());

          // Writing the decoded image as .dat file.
          if (result.count("output")) {
            matwrite(result["output"].as<std::string>() +
                         iterator->path().filename().string() + ".dat",
                     image, image_size_x, image_size_y);
          } else {
            matwrite(iterator->path().string() + ".dat", image, image_size_x,
                     image_size_y);
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
    unsigned int image_size_x = 0, image_size_y = 0;

    image = (int*)decoder.DecodeFile(result["file"].as<std::string>(),
                                     &image_size_x, &image_size_y,
                                     result["level"].as<int>());

    // Writing the decoded image as .dat file.
    matwrite(result["file"].as<std::string>() + ".dat", image, image_size_x,
             image_size_y);

    return 0;
  } else {
    std::cout << "No file specified, exiting." << std::endl;
  }

  return 0;
}
