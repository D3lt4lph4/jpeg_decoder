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

  this->frame_header_ = ParseFrameHeader(
      this->current_file_content_, this->current_index_, encoding_process_type);

  if (this->frame_header_.encoding_process_type_ == FRAME_TYPE_BASELINE_DTC) {
    this->number_of_blocks_per_column =
        (this->frame_header_.number_of_lines_ + 8 - 1) / 8;
    this->number_of_blocks_per_line =
        (this->frame_header_.number_of_samples_per_line_ + 8 - 1) / 8;
    this->current_image_ = cv::Mat(
        this->frame_header_.number_of_lines_ +
            (32 + (42 - (this->frame_header_.number_of_lines_ % 8)) % 8),
        this->frame_header_.number_of_samples_per_line_ +
            (42 - ((this->frame_header_.number_of_samples_per_line_ % 8)) % 8),
        CV_32SC3);
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
  this->scan_header_ =
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
  unsigned char vertical_number_of_blocks, dc_table_index, ac_table_index,
      bit_index = 8;
  unsigned int mcu_number, number_of_mcus, h_max = 0, v_max = 0;
  int prev[3] = {0, 0, 0};
  this->ResetDecoderBaseline();

  // Calculating the number of mcus to parse.

  for (size_t component_number = 1;
       component_number <= this->frame_header_.number_of_component_;
       component_number++) {
    if (this->frame_header_.component_parameters_.at(component_number).at(0) >
        h_max) {
      h_max =
          this->frame_header_.component_parameters_.at(component_number).at(0);
    }
    if (this->frame_header_.component_parameters_.at(component_number).at(1) >
        v_max) {
      v_max =
          this->frame_header_.component_parameters_.at(component_number).at(1);
    }
  }

  number_of_mcus =
      ((this->frame_header_.number_of_lines_ + (v_max * 8) - 1) / (v_max * 8)) *
      ((this->frame_header_.number_of_samples_per_line_ + (h_max * 8) - 1) /
       (h_max * 8));

  for (unsigned int mcu_number = 0; mcu_number < number_of_mcus; mcu_number++) {
    this->DecodeMCUBaseline(mcu_number, h_max, v_max, &bit_index, prev);
  }

  if (bit_index < 8) {
    *(this->current_index_) += 1;
  }
}

/**
 * \fn void JPEGDecoder::DecodeMCUBaseline(unsigned int mcu_number)
 * \brief decode a full mcu for the baseline encoding.
 * The function handles the level processing for each block of data processed.
 *
 * \param[in] mcu_number The index of the mcu being processed. The mcu indexing
 * starts from the top left mcu and goes line by line to the bottom right one.
 *
 */
void JPEGDecoder::DecodeMCUBaseline(unsigned int mcu_number, unsigned int h_max,
                                    unsigned int v_max,
                                    unsigned char *bit_index, int *prev) {
  unsigned int horizontal_number_of_blocks, vertical_number_of_blocks,
      start_line, start_column, number_of_mcu_per_line,
      number_of_mcu_per_column, mcu_per_line;

  int diff;

  unsigned char decoded_dc, dc_table_index, ac_table_index, number_of_component;

  cv::Mat new_block;
  std::vector<int> AC_Coefficients;

  number_of_component = this->frame_header_.number_of_component_;

  mcu_per_line =
      (this->frame_header_.number_of_samples_per_line_ + (h_max * 8) - 1) /
      (h_max * 8);
  start_line = mcu_number / mcu_per_line * v_max;
  start_column = mcu_number % mcu_per_line * h_max;

  for (unsigned char component_number = 1;
       component_number <= number_of_component; component_number++) {
    dc_table_index =
        this->scan_header_.components_parameters_.at(component_number).first;
    ac_table_index =
        this->scan_header_.components_parameters_.at(component_number).second;
    horizontal_number_of_blocks =
        this->frame_header_.component_parameters_.at(component_number).at(0);
    vertical_number_of_blocks =
        this->frame_header_.component_parameters_.at(component_number).at(1);

    // We process all the blocks for the current component.
    for (size_t v_block = 0; v_block < vertical_number_of_blocks; v_block++) {
      for (size_t h_block = 0; h_block < horizontal_number_of_blocks;
           h_block++) {
        // Getting the dc coefficient.
        decoded_dc =
            Decode(this->current_file_content_, this->current_index_, bit_index,
                   this->dc_huffman_tables_.at(dc_table_index));
        diff = Receive(decoded_dc, this->current_file_content_,
                       this->current_index_, bit_index);

        diff = Extended(diff, decoded_dc);
        prev[component_number - 1] = diff + prev[component_number - 1];

        // Getting the AC coefficients.
        AC_Coefficients = DecodeACCoefficients(
            this->current_file_content_, this->current_index_, bit_index,
            this->ac_huffman_tables_.at(ac_table_index));

        // We save the info in the correct blocks.
        for (size_t row_d = 0; row_d < v_max / vertical_number_of_blocks;
             row_d++) {
          for (size_t column_d = 0;
               column_d < h_max / horizontal_number_of_blocks; column_d++) {
            new_block = this->current_image_(
                cv::Range(8 * (start_line + row_d +
                               v_block * v_max / vertical_number_of_blocks),
                          8 * (1 + start_line + row_d +
                               v_block * v_max / vertical_number_of_blocks)),
                cv::Range(8 * (start_column + column_d +
                               h_block * h_max / horizontal_number_of_blocks),
                          8 * (1 + start_column + column_d +
                               h_block * h_max / horizontal_number_of_blocks)));

            // We save the dc coefficient and update the previous value.
            new_block.at<cv::Vec3i>(0, 0)[component_number - 1] =
                prev[component_number - 1];

            for (size_t row = 0; row < 8; row++) {
              for (size_t col = 0; col < 8; col++) {
                if (!(row == 0 && col == 0)) {
                  new_block.at<cv::Vec3i>(row, col)[component_number - 1] =
                      AC_Coefficients.at(ZZ_order[row * 8 + col] - 1);
                }
              }
            }
            // If required, dequantize the coefficient.
            if (this->decoding_level_ > 1) {
              // Perform dequantization
              this->Dequantize(&new_block,
                               this->quantization_tables_.at(
                                   this->frame_header_.component_parameters_
                                       .at((unsigned char)component_number)
                                       .at(2)),
                               component_number);
            }

            // If required Perform the dct inverse.
            if (this->decoding_level_ > 2) {
              IDCT(&new_block, component_number);
            }
          }
        }
      }
    }
    // If required, perform the color space tranform on the newly decoded
    // blocks.
  }
  if (this->decoding_level_ > 3) {
    for (size_t row_d = 0; row_d < v_max; row_d++) {
      for (size_t column_d = 0; column_d < h_max; column_d++) {
        new_block = this->current_image_(
            cv::Range(8 * (start_line + row_d), 8 * (1 + start_line + row_d)),
            cv::Range(8 * (start_column + column_d),
                      8 * (1 + start_column + column_d)));
        YCbCrToBGR(&new_block);
      }
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
 * The 00 are used to differanciate between the marker and values in the
 * encoded stream.
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
  int value;
  for (size_t i = 0; i < 8; i++) {
    for (size_t j = 0; j < 8; j++) {
      value = new_block->at<cv::Vec3i>(i, j)[component_number - 1];
      new_block->at<cv::Vec3i>(i, j)[component_number - 1] *=
          table.qks_.at(ZZ_order[i * 8 + j]);
    }
  }
}
