#include <iostream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGDecoder.hpp"

int main(int argc, char const *argv[])
{
    JPEGDecoder decoder;
    cv::Mat image;

    image = decoder.Decode("data/chat.jpg", 1);

    std::cout << decoder;

    return 0;
}
