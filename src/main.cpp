#include <iostream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGDecoder.hpp"

int main(int argc, char const *argv[])
{
    if (argc == 2)
    {
        JPEGDecoder decoder;
        cv::Mat image;

        image = decoder.Decode(argv[1], 1);

        cv::namedWindow("Decoded image.", cv::WINDOW_AUTOSIZE);
        cv::imshow("Decoded image.", image);
        cv::waitKey(0);
    }
    else
    {
        std::cout << "Please enter one image to decode." << std::endl;
    }

    return 0;
}
