#include <boost/filesystem.hpp>
#include <cerrno>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <cxxopts.hpp>

#include "JPEGDecoder.hpp"
#include "JPEGUtils.hpp"

// We only use opencv if we are in debug mode to visually check the results of
// the decoding
#ifdef DEBUG
#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#endif

void matwrite(const std::string &filename, JPEGImage *mat, int channels)
{
  std::ofstream fs(filename, std::fstream::binary);
  int number_of_component, temp;
  if (fs.fail())
  {
    std::cerr << strerror(errno) << std::endl;
  }

  // Writing the number of components
  number_of_component = mat->GetNumberOfComponent();
  fs.write((char *)&number_of_component, sizeof(int));

  // First we write all of the component sizes.
  for (size_t i = 0; i < number_of_component; i++)
  {
    temp = mat->GetComponentShape(i).first;
    fs.write((char *)&temp, sizeof(int));
    temp = mat->GetComponentShape(i).second;
    fs.write((char *)&temp, sizeof(int));
  }

  // We write the real image size, for future resizing if needed.
  for (size_t i = 0; i < number_of_component; i++)
  {
    temp = mat->GetRealShape().at(i);
    fs.write((char *)&temp, sizeof(int));
  }

  // Then we write all of the data.

  for (size_t i = 0; i < number_of_component; i++)
  {
    fs.write((char *)mat->GetData(i).data(),
             sizeof(int) * mat->GetComponentShape(i).first *
                 mat->GetComponentShape(i).second);
  }
  fs.close();
}

int main(int argc, char *argv[])
{
#ifdef DEBUG
  bool show = false;
#endif
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
  if (result.count("help"))
  {
    std::cout << options.help({"", "Group"}) << std::endl;
    exit(0);
  }

#ifdef DEBUG
  if (result.count("directory") && show)
  {
    std::cout << "Can only display for one image, exiting..." << std::endl;
    exit(0);
  }
#endif

  // If the directory is specified, we process the files in it.
  if (result.count("directory"))
  {
    std::string directory = result["directory"].as<std::string>();
    boost::filesystem::directory_iterator end_iterator;
    if (!boost::filesystem::exists(directory))
    {
      std::cout << "The directory doesn't exists, exiting..." << std::endl;
      exit(0);
    }

    for (boost::filesystem::directory_iterator iterator(directory);
         iterator != end_iterator; ++iterator)
    {
      if (!boost::filesystem::is_directory(iterator->status()))
      {
        if (boost::filesystem::extension(iterator->path()) == ".jpg")
        {
          JPEGDecoder decoder;
          JPEGImage *image;

          std::cout << "Processing the image : " << iterator->path().string()
                    << std::endl;
          image = decoder.DecodeFile(iterator->path().string(),
                                     result["level"].as<int>());

          // Writing the decoded image as .dat file.
          if (result.count("output"))
          {
            switch (result["level"].as<int>())
            {
            case 2:
              matwrite(result["output"].as<std::string>() +
                           iterator->path().stem().string() + ".qhjpg",
                       image, channels);
              break;
            case 3:
              matwrite(result["output"].as<std::string>() +
                           iterator->path().stem().string() + ".iqhjpg",
                       image, channels);
              break;
            default:
              matwrite(result["output"].as<std::string>() +
                           iterator->path().stem().string() + ".riqhjpg",
                       image, channels);
              break;
            }
          }
          else
          {
            switch (result["level"].as<int>())
            {
            case 2:
              matwrite(iterator->path().parent_path().string() + "/" +
                           iterator->path().stem().string() + ".qhjpg",
                       image, channels);
              break;
            case 3:
              matwrite(iterator->path().parent_path().string() + "/" +
                           iterator->path().stem().string() + ".iqhjpg",
                       image, channels);
              break;
            default:
              matwrite(iterator->path().parent_path().string() + "/" +
                           iterator->path().stem().string() + ".riqhjpg",
                       image, channels);
              break;
            }
          }
        }
      }
    }
    return 0;
  }

  // Process the file if specified.
  if (result.count("file"))
  {
    JPEGDecoder decoder(0);
    JPEGImage *image;
#ifdef DEBUG
    if (show)
    {
      std::cout << "Show is set, processing to level 4, full extraction."
                << std::endl;
      image = decoder.DecodeFile(result["file"].as<std::string>(), 4);
    }
    else
    {
      image = decoder.DecodeFile(result["file"].as<std::string>(),
                                 result["level"].as<int>());
    }
#else
    image = decoder.DecodeFile(result["file"].as<std::string>(),
                               result["level"].as<int>());
#endif

    // Writing the decoded image as .dat file.
    boost::filesystem::path p(result["file"].as<std::string>());

    if (result.count("output"))
    {
      switch (result["level"].as<int>())
      {
      case 2:
        matwrite(
            result["output"].as<std::string>() + p.stem().string() + ".qhjpg",
            image, channels);
        break;
      case 3:
        matwrite(result["output"].as<std::string>() + p.stem().string() +
                     ".iqhjpg",
                 image, channels);
        break;
      default:
        matwrite(result["output"].as<std::string>() + p.stem().string() +
                     ".riqhjpg",
                 image, channels);
        break;
      }
    }
    else
    {
      switch (result["level"].as<int>())
      {
      case 2:
        matwrite(
            p.parent_path().string() + "/" + p.stem().string() + ".qhjpg",
            image, channels);
        break;
      case 3:
        matwrite(
            p.parent_path().string() + "/" + p.stem().string() + ".iqhjpg",
            image, channels);
        break;
      default:
        matwrite(
            p.parent_path().string() + "/" + p.stem().string() + ".riqhjpg",
            image, channels);
        break;
      }
    }
#ifdef DEBUG
    if (show)
    {
      std::vector<int> shape = image->GetRealShape();
      cv::Mat *image_to_display;
      if (shape[2] == 1)
      {
        image_to_display = new cv::Mat(shape.at(0), shape.at(1), CV_32SC1);
      }
      else
      {
        image_to_display = new cv::Mat(shape.at(0), shape.at(1), CV_32SC3);
      }

      for (size_t row = 0; row < shape.at(0); row++)
      {
        for (size_t col = 0; col < shape.at(1); col++)
        {
          for (size_t channel = 0; channel < shape[2]; channel++)
          {
            if (shape[2] == 1)
            {
              image_to_display->at<int>(row, col) = image->at(row, col, channel);
            }
            else
            {
              image_to_display->at<cv::Vec3i>(row, col)[channel] = image->at(row, col, channel);
            }
          }
        }
      }

      image_to_display->convertTo(*image_to_display, CV_8UC3);
      cv::imshow("Decoded image.", *image_to_display);
      cv::waitKey(0);
    }
#endif
    return 0;
  }
  else
  {
    std::cout << "No file specified, exiting." << std::endl;
  }

  return 0;
}
