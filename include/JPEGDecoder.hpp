#ifndef __JPEG_DECODER__
#define __JPEG_DECODER__

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


class JPEGDecoder {

    public:
        cv::Mat Decode(std::string file_name, int level);


    private:

        void DecodeHuffman();
        void DecodeRunLengthEncoding();
        void InverseQuantification();
        void InverseDirectCosineTransform();

        bool GetFileInformation(std::ifstream file_to_decode);


        cv::Mat current_image_, current_thumbnail_;
        std::string current_filename_;

        
}
#endif