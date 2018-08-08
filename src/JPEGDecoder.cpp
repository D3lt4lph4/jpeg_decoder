#include <math.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"
#include "JPEGHuffmanDecoder.hpp"
#include "JPEGParser.hpp"
#include "JPEGUtility.hpp"

/**
 * \fn JPEGDecoder::JPEGDecoder()
 */
JPEGDecoder::JPEGDecoder() {
  this->current_index_ = new (unsigned int);
  this->block_index = 0;
}

/**
 * \fn cv::Mat JPEGDecoder::Decode(std::string filename, int level)
 * \brief Take a JPEG compress file as entry and output the decoded matrix.
 *
 * \param[in] filename The file to be decoded.
 * \param[in] level The required level of decoding.
 */
cv::Mat JPEGDecoder::DecodeFile(std::string filename, int level) {
  std::ifstream file_to_decode;
  int size, current_index = 0;
  unsigned char *marker, table_key;
  QuantizationTable quantization_table;
  HuffmanTable huffman_table;
  std::vector<std::pair<unsigned char, HuffmanTable>> huffman_tables;

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
      do {
        marker = this->GetMarker();
        switch (*marker) {
          case APPO:
            this->current_jfif_header = ParseJFIFSegment(
                this->current_file_content_, this->current_index_);
            std::cout << "JFIF segment parsed." << std::endl;
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
            std::tie(table_key, quantization_table) = ParseQuantizationTable(
                this->current_file_content_, this->current_index_);

            if (!(this->quantization_tables_
                      .insert(std::make_pair(table_key, quantization_table))
                      .second)) {
              this->quantization_tables_.erase(table_key);
              this->quantization_tables_.insert(
                  std::make_pair(table_key, quantization_table));
            }
            std::cout << "Quantization table parsed." << std::endl;
            break;
          case START_OF_FRAME_BASELINE:
            this->DecodeFrame(FRAME_TYPE_BASELINE_DTC);
            std::cout << "Frame decoded." << std::endl;
            break;
          case START_OF_FRAME_PROGRESSIVE:
            this->DecodeFrame(FRAME_TYPE_PROGRESSIVE);
            std::cout << "Frame decoded." << std::endl;
            break;
          case DEFINE_HUFFMAN_TABLE:
            huffman_tables = ParseHuffmanTableSpecification(
                this->current_file_content_, this->current_index_);

            for (size_t i = 0; i < huffman_tables.size(); i++) {
              std::tie(table_key, huffman_table) = huffman_tables.at(i);

              if (huffman_table.table_class_ == 0) {
                if (!(this->dc_huffman_tables_
                          .insert(std::make_pair(table_key, huffman_table))
                          .second)) {
                  this->dc_huffman_tables_.erase(table_key);
                  this->dc_huffman_tables_.insert(
                      std::make_pair(table_key, huffman_table));
                }
              } else {
                if (!(this->ac_huffman_tables_
                          .insert(std::make_pair(table_key, huffman_table))
                          .second)) {
                  this->ac_huffman_tables_.erase(table_key);
                  this->ac_huffman_tables_.insert(
                      std::make_pair(table_key, huffman_table));
                }
              }
            }
            std::cout << "Huffman table decoded." << std::endl;
            break;
          case END_OF_IMAGE:
            break;
          default:
            std::cout << "I did not know how to parse the block : " << std::hex
                      << int(*marker) << std::endl;
            throw std::runtime_error("Error while processing the jpeg file.");
            break;
        }
      } while (*marker != END_OF_IMAGE);
    } else {
      throw std::runtime_error(
          "Expected SOI marker, but sosmething else found, cannot parse the "
          "file.");
    }
  } else {
    throw FileNotFoundException(filename);
  }
  return this->current_image_;
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
  *(this->current_index_) = 0;
  this->block_index = 0;
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
  unsigned char *marker, table_key;
  QuantizationTable quantization_table;
  HuffmanTable huffman_table;
  std::vector<std::pair<unsigned char, HuffmanTable>> huffman_tables;
  ;

  this->current_frame_header_ = ParseFrameHeader(
      this->current_file_content_, this->current_index_, encoding_process_type);

  if (this->current_frame_header_.encoding_process_type_ ==
      FRAME_TYPE_BASELINE_DTC) {
    this->number_of_blocks_per_column =
        (this->current_frame_header_.number_of_lines_ + 8 - 1) / 8;
    this->number_of_blocks_per_line =
        (this->current_frame_header_.number_of_samples_per_line_ + 8 - 1) / 8;
    this->current_image_ = cv::Mat(
        this->current_frame_header_.number_of_lines_,
        this->current_frame_header_.number_of_samples_per_line_, CV_32FC3);
  }

  do {
    marker = this->GetMarker();
    if (*marker == START_OF_SCAN) {
      std::cout << "Getting scan" << std::endl;
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
          std::cout << "Quantization table parsed." << std::endl;
          break;
        case DEFINE_HUFFMAN_TABLE:
          huffman_tables = ParseHuffmanTableSpecification(
              this->current_file_content_, this->current_index_);
          for (size_t i = 0; i < huffman_tables.size(); i++) {
            std::tie(table_key, huffman_table) = huffman_tables.at(i);

            if (huffman_table.table_class_ == 0) {
              if (!(this->dc_huffman_tables_
                        .insert(std::make_pair(table_key, huffman_table))
                        .second)) {
                this->dc_huffman_tables_.erase(table_key);
                this->dc_huffman_tables_.insert(
                    std::make_pair(table_key, huffman_table));
              }
            } else {
              if (!(this->ac_huffman_tables_
                        .insert(std::make_pair(table_key, huffman_table))
                        .second)) {
                this->ac_huffman_tables_.erase(table_key);
                this->ac_huffman_tables_.insert(
                    std::make_pair(table_key, huffman_table));
              }
            }
          }
          std::cout << "Huffman table parsed." << std::endl;
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
  *(this->current_index_) -= 2;
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
  int n = 0, start_line, start_column;
  unsigned char decoded_dc;
  int diff;
  this->ResetDecoderBaseline();
  unsigned int component_number = 1, number_of_blocks;
  unsigned char dc_table_index, ac_table_index, bit_index = 8;
  this->data_unit_per_mcu_ = 3;
  cv::Mat new_block;
  std::vector<int> AC_Coefficients;
  number_of_blocks =
      this->number_of_blocks_per_line * this->number_of_blocks_per_column;
  while (this->block_index < number_of_blocks) {
    n = 0;
    while (n < this->data_unit_per_mcu_) {
      if (this->current_frame_header_.number_of_image_component > 1) {
        start_line = this->block_index / this->number_of_blocks_per_line;
        start_column = this->block_index % this->number_of_blocks_per_line;
        new_block = this->current_image_(
            cv::Range(8 * start_line, 8 * (1 + start_line)),
            cv::Range(8 * start_column, 8 * (1 + start_column)));
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
            Decode(this->current_file_content_, this->current_index_,
                   &bit_index, this->dc_huffman_tables_.at(dc_table_index));
        diff = Receive(decoded_dc, this->current_file_content_,
                       this->current_index_, &bit_index);

        diff = Extended(diff, decoded_dc);
        new_block.at<cv::Vec3i>(0, 0)[component_number - 1] = diff;
        // std::cout << diff << " : ";
        // We decode the ac components.
        AC_Coefficients = DecodeACCoefficients(
            this->current_file_content_, this->current_index_, &bit_index,
            this->ac_huffman_tables_.at(ac_table_index));

        if (this->decoding_level_ > 1) {
          // Perform dequantization
          this->Dequantize(&new_block,
                           this->quantization_tables_.at(
                               this->current_frame_header_
                                   .component_signification_parameters_
                                   .at((unsigned char)component_number)
                                   .at(2)),
                           component_number);
        }

        for (size_t i = 0; i < 8; i++) {
          for (size_t j = 0; j < 8; j++) {
            if (!(i == 0 && j == 0)) {
              new_block.at<cv::Vec3i>(i, j)[component_number - 1] =
                  AC_Coefficients.at(ZZ_order[i * 8 + j - 1]);
              // std::cout << AC_Coefficients.at(ZZ_order[i * 8 + j - 1]);
              // std::cout << " ";
            }
          }
          // std::cout << std::endl;
        }

        if (this->decoding_level_ > 2) {
          // Perform IDCT.
          IDCT(&new_block, component_number);
        }

        if (component_number ==
            this->current_frame_header_.number_of_image_component) {
          component_number = 1;
          this->block_index += 1;
        } else {
          component_number += 1;
        }
      }
      if (this->decoding_level_ > 3) {
        // Perform IDCT.
        YCbCrToBGR(&new_block);
      }
      n += 1;
    }
  }
  if (bit_index < 8) {
    *(this->current_index_) += 1;
  }
}

/**
 * \fn bool JPEGDecoder::IsMarker()
 * \brief Tells whether the next two bits in the stream are of a marker.
 * This function does not increment the index count for the stream processed.
 */
bool JPEGDecoder::IsMarker() {
  if (this->current_file_content_[*(this->current_index_)] == 0xFF) {
    if (this->current_file_content_[*(this->current_index_) + 1] == 0x00) {
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
    *(this->current_index_) = *(this->current_index_) + 2;
    return marker;
  } else {
    std::cout << this->current_index_ << std::endl;
    error << "Error while reading marker, 0xFF expected, but " << std::hex
          << (int)this->current_file_content_[*(this->current_index_)]
          << " found.";
    throw std::runtime_error(error.str());
  }
}

void JPEGDecoder::Dequantize(cv::Mat *new_block, QuantizationTable table,
                             unsigned int component_number) {
  for (size_t i = 0; i < 8; i++) {
    for (size_t j = 0; j < 8; j++) {
      new_block->at<cv::Vec3i>(i, j)[component_number - 1] *=
          table.qks_.at(i * 8 + j);
    }
  }
}
