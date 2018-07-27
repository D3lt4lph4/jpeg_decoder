#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"

cv::Mat JPEGDecoder::Decode(std::string filename, int level) {

    std::ifstream file_to_decode;
    int size, current_index = 0;
    unsigned char* file_content, *marker;
    // If the filename is the same, we assume to have the same image.
    if (filename.compare(this->current_filename_) == 0) {
        return this->current_image_.clone();
    }

    file_to_decode.open(filename, std::ios::binary);

    if (file_to_decode.is_open()) {
        file_to_decode.seekg(0, std::ios::end);
        size = file_to_decode.tellg();
        file_to_decode.seekg(0, std::ios::beg);
        file_content = new unsigned char[size + 1];
        file_to_decode.read(reinterpret_cast<char*>(file_content), size);

        if (this->GetFileInformation(file_content, &current_index)) {
            while (file_content[current_index] != END_OF_IMAGE) {
                marker = this->GetMarker(file_content,  &current_index, 2);
                
                switch (*marker)
                {
                    case DEFINE_QUANTIZATION_TABLE:
                        this->ParseQuantizationTable(file_content, &current_index);
                        break;
                
                    default:
                        throw std::runtime_error("Error while processing the jpeg file.");
                        break;
                }
            }

        } else {
            throw std::runtime_error("Error while getting the information on the file to parse.");
        }
    } else {
        throw FileNotFoundException(filename);
    }

    return cv::Mat(300, 300, 25);

}

bool JPEGDecoder::GetFileInformation(unsigned char *file_content, int *index) {

    unsigned char* marker;

    // Start of image marker.
    marker = this->GetMarker(file_content, index, 2);
    if (*marker != START_OF_IMAGE) {
        return false;
    }

    // Getting APP0 marker.
    marker = this->GetMarker(file_content, index, 2);

    if (*marker != APPO) {
        
        return false;
    }

    // Getting JFIF marker
    marker = this->GetMarker(file_content, index, 2, false);
    marker = this->GetMarker(file_content, index, 5, false);
    if (memcmp(marker, JFIF, 5) != 0) {
        return false;
    }
    
    // Getting the version of the encoding.
    marker = this->GetMarker(file_content, index, 2, false);
    this->current_version_ = int(marker[0] << 8 | marker[1]);

    // Getting the units.
    marker = this->GetMarker(file_content, index, 1, false);
    this->current_unit_ = int(marker[0]);

    // Getting the Horizontal Pixel Density
    marker = this->GetMarker(file_content, index, 2, false);
    this->horizontal_pixel_density_ = int(marker[0] << 8 | marker[1]);

    // Getting the Vertical Pixel Density
    marker = this->GetMarker(file_content, index, 2, false);
    this->vertical_pixel_density_ = int(marker[0] << 8 | marker[1]);

    // Getting the Thumbnail Horizontal Pixel Count
    marker = this->GetMarker(file_content, index, 1, false);
    this->thumbnail_horizontal_pixel_count_ = int(marker[0]);

    // Getting the Thumbnail Vertical Pixel Count
    marker = this->GetMarker(file_content, index, 1, false);
    this->thumbnail_vertical_pixel_count_ = int(marker[0]);

    // If we have thumbnail count, we get the image.
    if (this->thumbnail_horizontal_pixel_count_ != 0 && this->thumbnail_vertical_pixel_count_ != 0) {
        std::cout << "TODO" << std::endl;
    }

    return true;
}

unsigned char* JPEGDecoder::GetMarker(unsigned char *file_content, int *index, int size, bool FF_expected) {
    std::stringstream error;
    unsigned char* marker;
    
    if (FF_expected) {
        marker = new unsigned char [size -1];
        std::memcpy(marker, &file_content[*index + 1], size - 1);
        
        if (file_content[*index] == 0xFF) {
            *index = *index + 2;
            return marker;
        } else {
            error << "Error while reading marker, 0xFF expected, but " << std::hex << (int) file_content[*index] << " found.";
            throw std::runtime_error(error.str());
        }
    }
    marker = new unsigned char [size];
    std::memcpy(marker, &file_content[*index], size);
    *index += size;
    return marker;
}

std::ostream& operator<<(std::ostream &out, const JPEGDecoder &decoder) {
    out << "Format: JFIF\n";
    out << "Version : " << decoder.current_version_ << "\n";
    
    switch (decoder.current_unit_)
    {
        case 0:
            out << "Unit : 0, no units\n";
            break;
        
        case 1:
            out << "Unit : 1, dot per inch\n";
            break;
        
        case 2:
            out << "Unit : 2, dot per cm\n";
            break;
        default:
            out << "No unit selected, this should not be possible.\n";
            break;
    }

    out << "Horizontal pixel density : " << decoder.horizontal_pixel_density_ << "\n";
    out << "Vertical pixel density : " << decoder.vertical_pixel_density_ << "\n";
    out << "Thumbnail horizontal pixel count : " << decoder.thumbnail_horizontal_pixel_count_ << "\n";
    out << "Thumbnail vertical pixel count : " << decoder.thumbnail_vertical_pixel_count_ << "\n";
    return out;
}

JPEGDecoder::ParseQuantizationTable(unsigned char)