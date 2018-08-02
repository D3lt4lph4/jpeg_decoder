#include <iostream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGDecoder.hpp"

int main(int argc, char const *argv[])
{
    JPEGDecoder decoder;
    cv::Mat image;

    image = decoder.Decode("data/private.jpg", 1);

    cv::namedWindow("Decoded image.", cv::WINDOW_AUTOSIZE);
    cv::imshow("Decoded image.", image);
    cv::waitKey(0);

    return 0;
}
