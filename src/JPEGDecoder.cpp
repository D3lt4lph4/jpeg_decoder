#include <string>
#include <stdexcept>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"

JPEGDecoder::Decode(std::string filename, int level) {

    std::ifstream file_to_decode;

    // If the filename is the same, we assume to have the same image.
    if (filename.compare(this.current_filename_) == 0) {
        return this.current_image.clone()
    }

    this.file_to_decode.open(filename, std::ios::binary)

    if (file_to_decode.is_open()) {
        if (this.GetFileInformation(file_to_decode)) {

        } else {
            throw runtime_error("Error while getting the informaiton on the file to parse.")
        }
    } else {
        throw FileNotFoundException("Could not open file to parse.");
    }

}

JPEGDecoder::GetFileInformation(std::ifstream file_to_decode) {
    unsigned char high, low;

    file_to_decode.read(reinterpret_cast<char*>(&high))
    file_to_decode.read(reinterpret_cast<char*>(&low))

    if (high != 0xFF || low != 0xD8) {
        return false;
    }

    file_to_decode.read(reinterpret_cast<char*>(&high))
    file_to_decode.read(reinterpret_cast<char*>(&low))

    if (high != 0xFF || low != 0xE0) {
        return false;
    }
    return true;
}