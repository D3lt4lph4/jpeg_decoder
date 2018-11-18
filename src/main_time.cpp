#include <boost/filesystem.hpp>
#include <cerrno>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

#include <cxxopts.hpp>

#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>

#include "JPEGDecoder.hpp"
#include "JPEGUtils.hpp"

int main(int argc, char* argv[]) {
  int channels;

  // Creating the parser.
  cxxopts::Options options("Time calculator",
                           "Program to calculate the difference in calculation "
                           "between OpenCV and the jpegdecoder.");
  options.add_options()(
      "d,directory",
      "The directory containing all the images to use for the test.",
      cxxopts::value<std::string>()->default_value("default"))("help",
                                                               "Print help");

  auto result = options.parse(argc, argv);

  // If the help is required, display it.
  if (result.count("help")) {
    std::cout << options.help({"", "Group"}) << std::endl;
    exit(0);
  }
  
  // If we have the argument (else we return a warning message).
  if (result.count("directory")) {
    std::string directory = result["directory"].as<std::string>();
    boost::filesystem::directory_iterator end_iterator;
    if (!boost::filesystem::exists(directory)) {
      std::cout << "The directory doesn't exists, exiting..." << std::endl;
      exit(0);
    }

    // We start the processing of the images with OpenCV (RGB)
    std::cout << "Starting to process the images with OpenCV." << std::endl;

    const clock_t begin_time = std::clock();
    for (boost::filesystem::directory_iterator iterator(directory);
         iterator != end_iterator; ++iterator) {
      if (!boost::filesystem::is_directory(iterator->status())) {
        if (boost::filesystem::extension(iterator->path()) == ".jpg") {
          JPEGDecoder decoder;
          JPEGImage* image;

          std::cout << "Processing the image : " << iterator->path().string()
                    << std::endl;
          image = decoder.DecodeFile(iterator->path().string(),
                                     result["level"].as<int>());

        }
      }
    }
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

    // We start the processing of the images with the jpegdecoder (RGB)
    std::cout << "Starting to process the images with the jpegdecoder (RGB)." << std::endl;
    const clock_t begin_time = std::clock();  
    for (boost::filesystem::directory_iterator iterator(directory);
         iterator != end_iterator; ++iterator) {
      if (!boost::filesystem::is_directory(iterator->status())) {
        if (boost::filesystem::extension(iterator->path()) == ".jpg") {
          JPEGDecoder decoder;
          JPEGImage* image;

          std::cout << "Processing the image : " << iterator->path().string()
                    << std::endl;
          image = decoder.DecodeFile(iterator->path().string(),
                                     result["level"].as<int>());

        }
      }
    }
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    
    // We start the processing of the images with the jpegdecoder (DCT)
    std::cout << "Starting to process the images with the jpegdecoder (DCT)." << std::endl;
    const clock_t begin_time = std::clock();  
    for (boost::filesystem::directory_iterator iterator(directory);
         iterator != end_iterator; ++iterator) {
      if (!boost::filesystem::is_directory(iterator->status())) {
        if (boost::filesystem::extension(iterator->path()) == ".jpg") {
          JPEGDecoder decoder;
          JPEGImage* image;

          std::cout << "Processing the image : " << iterator->path().string()
                    << std::endl;
          image = decoder.DecodeFile(iterator->path().string(),
                                     result["level"].as<int>());

        }
      }
    }
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

  } else {
    std::cout << "You must specify a directory for me to calculate some time." << std::cout;
  }


  return 0;
}
