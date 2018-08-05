#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"
#include "JPEGHuffmanDecoder.hpp"
#include "JPEGParser.hpp"

/**
 * \fn JPEGDecoder::JPEGDecoder()
 */
JPEGDecoder::JPEGDecoder() {}

/**
 * \fn cv::Mat JPEGDecoder::Decode(std::string filename, int level)
 * \brief Take a JPEG compress file as entry and output the decoded matrix.
 *
 * \param[in] filename The file to be decoded.
 * \param[in] level The required level of decoding.
 */
cv::Mat JPEGDecoder::Decode(std::string filename, int level) {
  std::ifstream file_to_decode;
  int size, current_index = 0;
  unsigned char *marker;

  // If the filename is the same, we assume to have the same image.
  if (filename.compare(this->current_filename_) == 0) {
    return this->current_image_.clone();
  }

  this->InitializeDecoder();
  this->decoding_level_ = level;

  file_to_decode.open(filename, std::ios::binary);

  if (file_to_decode.is_open()) {
    // Getting the file content to be read.
    file_to_decode.seekg(0, std::ios::end);
    size = file_to_decode.tellg();
    file_to_decode.seekg(0, std::ios::beg);
    this->current_file_content_ = new unsigned char[size + 1];
    file_to_decode.read(reinterpret_cast<char *>(this->current_file_content_),
                        size);

    if (*(this->GetMarker()) == START_OF_IMAGE) {
      this->DecoderSetup();
      while (this->current_file_content_[current_index] != END_OF_IMAGE) {
        marker = this->GetMarker();

        switch (*marker) {
          case APPO:
            this->current_jfif_header = ParseJFIFSegment(
                this->current_file_content_, this->current_index_);
            break;
          case COMMENT:
            std::cout << ParseComment(this->current_file_content_,
                                      this->current_index_)
                      << std::endl;
            break;
          case DEFINE_RESTART_INTERVAL:
            // TODO
            break;
          case DEFINE_QUANTIZATION_TABLE:
            ParseQuantizationTable(this->current_file_content_,
                                   this->current_index_);
            break;
          case START_OF_FRAME_BASELINE:
            this->DecodeFrame(FRAME_TYPE_BASELINE_DTC);
            break;
          case START_OF_FRAME_PROGRESSIVE:
            this->DecodeFrame(FRAME_TYPE_PROGRESSIVE);
            break;
          case DEFINE_HUFFMAN_TABLE:
            ParseHuffmanTableSpecification(this->current_file_content_,
                                           this->current_index_);
            break;
          default:
            std::cout << "I did not know how to parse the block : " << std::hex
                      << int(*marker) << std::endl;
            throw std::runtime_error("Error while processing the jpeg file.");
            break;
        }
      }
    } else {
      throw std::runtime_error(
          "Expected SOI marker, but sosmething else found, cannot parse the "
          "file.");
    }
  } else {
    throw FileNotFoundException(filename);
  }

  return cv::Mat(300, 300, 25);
}

std::ostream &operator<<(std::ostream &out, const JPEGDecoder &decoder) {
  out << "Format: JFIF\n";
  out << "Version : " << decoder.current_jfif_header.current_version_ << "\n";

  switch (decoder.current_jfif_header.current_unit_) {
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

  out << "Horizontal pixel density : "
      << decoder.current_jfif_header.horizontal_pixel_density_ << "\n";
  out << "Vertical pixel density : "
      << decoder.current_jfif_header.vertical_pixel_density_ << "\n";
  out << "Thumbnail horizontal pixel count : "
      << decoder.current_jfif_header.thumbnail_horizontal_pixel_count_ << "\n";
  out << "Thumbnail vertical pixel count : "
      << decoder.current_jfif_header.thumbnail_vertical_pixel_count_ << "\n";
  return out;
}

/**
 * \fn void JPEGDecoder::InitializeDecoder()
 * \brief Initialize all the values in the decoder for a new image to decode.
 */
void JPEGDecoder::InitializeDecoder() {
  // We remove all elements from the maps of tables if any.
  this->quantization_tables_.clear();
  this->dc_huffman_tables_.clear();
  this->ac_huffman_tables_.clear();
}

/**
 * \fn void JPEGDecoder::DecoderSetup()
 * \brief Sets up the decoder.
 */
void JPEGDecoder::DecoderSetup() { this->restart_interval = 0; }

/**
 * \fn void JPEGDecoder::DecodeFrame(unsigned char encoding_process_type)
 * \brief Decode a frame. Is called by the decode function when a frame marker
 * is spotted.
 *
 * \param[in] encoding_process_type The type of the frame to decode.
 */
void JPEGDecoder::DecodeFrame(unsigned char encoding_process_type) {
  unsigned char *marker;

  ParseFrameHeader(this->current_file_content_, this->current_index_,
                   encoding_process_type);

  do {
    marker = this->GetMarker();
    if (*marker == START_OF_SCAN) {
      this->DecodeScan(encoding_process_type);
    } else {
      switch (*marker) {
        case COMMENT:
          ParseComment(this->current_file_content_, this->current_index_);
          break;
        case DEFINE_RESTART_INTERVAL:
          // TODO
          break;
        case DEFINE_QUANTIZATION_TABLE:
          ParseQuantizationTable(this->current_file_content_,
                                 this->current_index_);
          break;
        case DEFINE_HUFFMAN_TABLE:
          ParseHuffmanTableSpecification(this->current_file_content_,
                                         this->current_index_);
          break;
        case END_OF_IMAGE:
          break;
        default:
          std::cout << "I did not know how to parse the block : " << std::hex
                    << int(*marker) << std::endl;
          throw std::runtime_error("Error while processing the jpeg file.");
          break;
      }
    }
  } while (*marker != END_OF_IMAGE);
}

/**
 * \fn void JPEGDecoder::DecodeScan(unsigned char encoding_process_type)
 * \brief Decode a scan. Is called by Decode when a scan marker is
 * spotted.
 *
 * \param[in] encoding_process_type The type of encoding currently used.
 */
void JPEGDecoder::DecodeScan(unsigned char encoding_process_type) {
  this->current_scan_ =
      ParseScanHeader(this->current_file_content_, this->current_index_);
  unsigned int m = 0;

  do {
    m += 1;
    switch (encoding_process_type) {
      case FRAME_TYPE_BASELINE_DTC:
        this->DecodeRestartIntervalBaseline();
        break;
      default:
        break;
    }
  } while (m < this->restart_interval);
}

/**
 * \fn void JPEGDecoder::ResetDecoderBaseline()
 */
void JPEGDecoder::ResetDecoderBaseline() {}

/**
 * \fn void JPEGDecoder::DecodeRestartIntervalBaseline()
 * \brief Parse coded blocks of information until the end of the restart
 * interval.
 *
 * For baseline with multiple components, the components are interleaved,
 * meaning the blocks will be coded in this order:
 *  - C1,1 ; C1,2 ; C1,3
 *
 * With C1,2 being the second component of the first block.
 *
 * Depending on the ratio between the component can be ordered a bit
 * differently, with more blocks or less blocks depending on the component.
 */
void JPEGDecoder::DecodeRestartIntervalBaseline() {
  int n = 0;
  unsigned char decoded_dc, diff;
  this->ResetDecoderBaseline();
  unsigned int component_number = 1;
  unsigned char dc_table_index, ac_table_index;
  this->data_unit_per_mcu_ = 3;
  cv::Mat new_block = cv::Mat::zeros(8, 8, CV_32FC1);

  while (!this->IsMarker()) {
    n = 0;
    while (n < this->data_unit_per_mcu_) {
      if (this->current_frame_header_.number_of_image_component > 1) {
        // Then we are interleaved.
        dc_table_index =
            this->current_scan_.scan_components_specification_parameters_
                .at((unsigned char)component_number)
                .first;
        ac_table_index =
            this->current_scan_.scan_components_specification_parameters_
                .at((unsigned char)component_number)
                .second;

        // We decode the DC component.
        decoded_dc =
            DecodeBaseline(this->dc_huffman_tables_.at(dc_table_index));
        diff = ReceiveBaseline(decoded_dc);
        diff = ExtendedBaseline(diff, decoded_dc);

        new_block.at<uchar>(0, 0, 0) = diff;

        // We decode the ac components.
        DecodeACCoefficients(&new_block,
                             this->ac_huffman_tables_.at(ac_table_index));

        if (component_number ==
            this->current_frame_header_.number_of_image_component) {
          component_number = 1;
        } else {
          component_number += 1;
        }
      }
      n += 1;
    }
  }
}

/**
 * \fn bool JPEGDecoder::IsMarker()
 * \brief Tells whether the next two bits in the stream are of a marker.
 * This function does not increment the index count for the stream processed.
 */
bool JPEGDecoder::IsMarker() {
  if (this->current_file_content_[*(this->current_index_)] == 0xFF) {
    if (this->current_file_content_[*(this->current_index_) + 1] == 0xFF) {
      return false;
    }
    return true;
  } else {
    return false;
  }
}

/**
 * \fn unsigned char *JPEGDecoder::GetMarker()
 * \brief Returns the current marker if any, else throws a runtime error.
 *
 * All of the marker are of the form 0xFF.. with the .. being anything but 00.
 * The 00 are used to differanciate between the marker and values in the encoded
 * stream.
 *
 * This function expect to get a marker, it will throw an error if none found.
 * This function increments the index count for the stream processed.
 */
unsigned char *JPEGDecoder::GetMarker() {
  std::stringstream error;
  unsigned char *marker;

  if (this->current_file_content_[*(this->current_index_)] == 0xFF) {
    marker = new unsigned char[1];
    std::memcpy(marker,
                &(this->current_file_content_[*(this->current_index_) + 1]), 1);
    this->current_index_ = this->current_index_ + 2;
    return marker;
  } else {
    std::cout << this->current_index_ << std::endl;
    error << "Error while reading marker, 0xFF expected, but " << std::hex
          << (int)this->current_file_content_[*(this->current_index_)]
          << " found.";
    throw std::runtime_error(error.str());
  }
}