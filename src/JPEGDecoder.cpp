/**
 * \file JPEGDecoder.cpp
 * The source code for the class JPEGDecoder
 *
 * \author Benjamin Deguerre
 * \version 1.0
 *
 */

#include <math.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"
#include "JPEGHuffmanDecoder.hpp"
#include "JPEGParser.hpp"
#include "JPEGUtility.hpp"

/** \class JPEGDecoder
 * \brief Decoding class
 */

/**
 * \fn JPEGDecoder::JPEGDecoder()
 */
JPEGDecoder::JPEGDecoder()
    : block_index(0), logging_level_(0), current_index_(0) {
  this->InitializeLogger();
  this->current_image_ = NULL;
}

/**
 * \fn JPEGDecoder::JPEGDecoder(unsigned char logging_level)
 *
 * \param[in] logging_level The level of logs to display.
 */
JPEGDecoder::JPEGDecoder(const unsigned char logging_level)
    : block_index(0), current_index_(0) {
  this->logging_level_ = logging_level;
  this->current_image_ = NULL;
}

JPEGDecoder::~JPEGDecoder() {
  if (this->current_image_ != NULL) {
    delete this->current_image_;
  }
}

/**
 * \fn JPEGImage *JPEGDecoder::DecodeFile(const std::string filename, const int
 * level)
 *
 * \brief Take a JPEG compress file as entry and output the decoded image at the
 * required level.
 *
 * \param[in] filename The file to be decoded.
 * \param[in] level The required level of decoding.
 */
JPEGImage *JPEGDecoder::DecodeFile(const std::string filename,
                                   const int level) {
  std::ifstream file_to_decode;
  int size, current_index = 0;
  unsigned char table_key;
  std::unique_ptr<unsigned char> marker;

  // When to stop the reading of the file
  bool out_condition = false;

  // The tables
  QuantizationTable quantization_table;
  HuffmanTable huffman_table;
  std::vector<std::pair<unsigned char, HuffmanTable>> huffman_tables;

  // If the filename is the same, we assume to have the same image.
  if (filename.compare(this->current_filename_) == 0) {
    return this->current_image_;
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
    this->current_index_ = 0;
    file_to_decode.read(reinterpret_cast<char *>(this->current_file_content_),
                        size);
    marker = this->GetMarker();
    if (*marker == START_OF_IMAGE) {
      this->DecoderSetup();
      do {
        marker = this->GetMarker();
        switch (*marker) {
          case APPO:
            this->current_jfif_header = ParseJFIFSegment(
                this->current_file_content_, this->current_index_);
            BOOST_LOG_TRIVIAL(info) << "JFIF segment parsed.";
            break;
          case APP1:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP2:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP3:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP4:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP5:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP6:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP7:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP8:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP9:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP10:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP11:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP12:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP13:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP14:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case APP15:
            BOOST_LOG_TRIVIAL(info) << "Application block 1 found, ignoring.";
            ParseApplicationBlock(this->current_file_content_,
                                  this->current_index_);
            break;
          case COMMENT: {
            std::string comment;
            comment =
                ParseComment(this->current_file_content_, this->current_index_);
            BOOST_LOG_TRIVIAL(info) << "Comment in the file : " << comment;
          } break;
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
            BOOST_LOG_TRIVIAL(info) << "Quantization table parsed.";
            break;
          case START_OF_FRAME_BASELINE:
            this->DecodeFrame(FRAME_TYPE_BASELINE_DTC);
            BOOST_LOG_TRIVIAL(info) << "Frame decoded.";
            break;
          case START_OF_FRAME_PROGRESSIVE:
            this->DecodeFrame(FRAME_TYPE_PROGRESSIVE);
            BOOST_LOG_TRIVIAL(info) << "Frame decoded.";
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
            BOOST_LOG_TRIVIAL(info) << "Huffman table decoded.";
            break;
          case END_OF_IMAGE:
            out_condition = true;
            break;
          default:
            BOOST_LOG_TRIVIAL(error)
                << "I did not know how to parse the block : " << std::hex
                << int(*marker);
            throw std::runtime_error("Error while processing the jpeg file.");
            break;
        }
      } while (!out_condition);
    } else {
      throw std::runtime_error(
          "Expected SOI marker, but sosmething else found, cannot parse the "
          "file.");
    }
  } else {
    throw FileNotFoundException(filename);
  }

  if (this->decoding_level_ > 3) {
    DeLevelShift(*(this->current_image_));
    this->current_image_->RescaleToRealSize();
    YCbCrToBGR(*(this->current_image_), this->current_image_->GetRealShape());
  }

  delete[] this->current_file_content_;
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
  this->current_index_ = 0;
  this->block_index = 0;
}

/**
 * \fn void JPEGDecoder::DecoderSetup()
 * \brief Sets up the decoder.
 */
void JPEGDecoder::DecoderSetup() { this->restart_interval = 0; }

/**
 * \fn void JPEGDecoder::DecodeFrame(const unsigned char encoding_process_type)
 * \brief Decode a frame. Is called by the decode function when a frame marker
 * is spotted.
 *
 * \param[in] encoding_process_type The type of the frame to decode. Only
 * baseline is supported for now.
 */
void JPEGDecoder::DecodeFrame(const unsigned char encoding_process_type) {
  unsigned char table_key;
  std::unique_ptr<unsigned char> marker;

  QuantizationTable quantization_table;
  HuffmanTable huffman_table;
  std::vector<std::pair<unsigned char, HuffmanTable>> huffman_tables;
  ;
  unsigned int h_max = 0, v_max = 0;

  this->frame_header_ = ParseFrameHeader(
      this->current_file_content_, this->current_index_, encoding_process_type);

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
  std::vector<std::pair<int, int>> sizes(
      this->frame_header_.number_of_component_);

  if (this->frame_header_.encoding_process_type_ == FRAME_TYPE_BASELINE_DTC) {
    this->number_of_blocks_per_column = (this->frame_header_.number_of_lines_);
    this->number_of_blocks_per_line =
        (this->frame_header_.number_of_samples_per_line_);

    std::vector<int> realShape(3);

    int size_factor_v =
        (((this->number_of_blocks_per_column + 8 - 1) / 8) + v_max - 1) /
        (v_max);
    int size_factor_h =
        (((this->number_of_blocks_per_line + 8 - 1) / 8) + h_max - 1) / (h_max);
    for (size_t component_number = 1;
         component_number <= this->frame_header_.number_of_component_;
         component_number++) {
      if (component_number == 1) {
        realShape[0] = this->number_of_blocks_per_column;
        realShape[1] = this->number_of_blocks_per_line;
        realShape[2] = 3;
      }
      sizes[component_number - 1].second =
          size_factor_h *
          this->frame_header_.component_parameters_.at(component_number).at(0) *
          8;
      sizes[component_number - 1].first =
          size_factor_v *
          this->frame_header_.component_parameters_.at(component_number).at(1) *
          8;
    }

    if (this->current_image_ != NULL) {
      delete this->current_image_;
    }

    this->current_image_ = new JPEGImage(sizes);
    this->current_image_->SetRealShape(realShape);
  }

  do {
    marker = this->GetMarker();
    if (*marker == START_OF_SCAN) {
      BOOST_LOG_TRIVIAL(info) << "Getting a scan to decode.";
      this->DecodeScan(encoding_process_type);
    } else {
      switch (*marker) {
        case COMMENT:
          BOOST_LOG_TRIVIAL(info) << "Comment in the current scan : "
                                  << ParseComment(this->current_file_content_,
                                                  this->current_index_);
          break;
        case DEFINE_RESTART_INTERVAL:
          // TODO
          break;
        case DEFINE_QUANTIZATION_TABLE:
          ParseQuantizationTable(this->current_file_content_,
                                 this->current_index_);
          BOOST_LOG_TRIVIAL(info) << "Quantization table parsed.";
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
          BOOST_LOG_TRIVIAL(info) << "Huffman table parsed.";
          break;
        case END_OF_IMAGE:
          break;
        default:
          BOOST_LOG_TRIVIAL(error)
              << "I did not know how to parse the block : " << std::hex
              << int(*marker);
          throw std::runtime_error("Error while processing the jpeg file.");
          break;
      }
    }

  } while (*marker != END_OF_IMAGE);
  this->current_index_ -= 2;
}

/**
 * \fn void JPEGDecoder::DecodeScan(const unsigned char encoding_process_type)
 * \brief Decode a scan. Is called by Decode when a scan marker is spotted.
 *
 * \param[in] encoding_process_type The type of encoding currently used.
 */
void JPEGDecoder::DecodeScan(const unsigned char encoding_process_type) {
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
  std::vector<int> prev = {0, 0, 0};
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
    this->DecodeMCUBaseline(mcu_number, h_max, v_max, bit_index, prev);
  }
  unsigned char marker[2];

  marker[0] = this->current_file_content_[this->current_index_];
  marker[1] = this->current_file_content_[this->current_index_ + 1];

  if (marker[0] == 0xFF && marker[1] == 0x00) {
    this->current_index_ += 2;
  } else if (bit_index < 8) {
    this->current_index_ += 1;
  }
}

/**
 * \fn void JPEGDecoder::DecodeMCUBaseline(const unsigned int mcu_number, const
 * unsigned int h_max, const unsigned int v_max, unsigned char &bit_index,
 * std::vector<int> &prev)
 *
 * \brief decode a full mcu for the baseline encoding. The function handles the
 * level processing for each block of data processed.
 *
 * \param[in] mcu_number The index of the mcu being processed. The mcu indexing
 * starts from the top left mcu and goes line by line to the bottom right one.
 *
 * \param[in] h_max The biggest value for h.
 *
 * \param[in] v_max The biggest value for v.
 *
 * \param[in, out] bit_index The index of the current bit in the byte being
 * parsed.
 *
 * \param[in, out] prev A vector holding the previous DC values for each
 * component.
 *
 */
void JPEGDecoder::DecodeMCUBaseline(const unsigned int mcu_number,
                                    const unsigned int h_max,
                                    const unsigned int v_max,
                                    unsigned char &bit_index,
                                    std::vector<int> &prev) {
  unsigned int horizontal_number_of_blocks, vertical_number_of_blocks,
      start_line, start_column, number_of_mcu_per_line,
      number_of_mcu_per_column, mcu_per_line;

  int diff;

  unsigned char decoded_dc, dc_table_index, ac_table_index, number_of_component;

  int line_length;
  std::vector<int> AC_Coefficients;

  number_of_component = this->frame_header_.number_of_component_;

  mcu_per_line =
      (this->frame_header_.number_of_samples_per_line_ + (h_max * 8) - 1) /
      (h_max * 8);

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

    line_length =
        mcu_per_line *
        this->frame_header_.component_parameters_.at(component_number).at(0) *
        8;
    start_line =
        mcu_number / mcu_per_line *
        this->frame_header_.component_parameters_.at(component_number).at(1) *
        8;
    start_column =
        mcu_number % mcu_per_line *
        this->frame_header_.component_parameters_.at(component_number).at(0) *
        8;
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

        // We save the dc coefficient and update the previous value.
        this->current_image_->at(
            start_line + 8 * v_block, start_column + 8 * h_block,
            component_number - 1) = prev[component_number - 1];

        for (size_t row = 0; row < 8; row++) {
          for (size_t col = 0; col < 8; col++) {
            if (!(row == 0 && col == 0)) {
#ifdef DEBUG
              this->current_image_->at(start_line + 8 * v_block + row,
                                       start_column + 8 * h_block + col,
                                       component_number - 1) =
                  AC_Coefficients.at(ZZ_order[row * 8 + col] - 1);
#else
              this->current_image_->at(start_line + 8 * v_block + row,
                                       start_column + 8 * h_block + col,
                                       component_number - 1) =
                  AC_Coefficients[ZZ_order[row * 8 + col] - 1];
#endif
            }
          }
        }
        // If required, dequantize the coefficient.
        if (this->decoding_level_ > 1) {
          // Perform dequantization
          this->Dequantize(component_number - 1, start_line + 8 * v_block,
                           start_column + 8 * h_block,
                           this->quantization_tables_.at(
                               this->frame_header_.component_parameters_
                                   .at((unsigned char)component_number)
                                   .at(2)));
        }

        // If required Perform the dct inverse.
        if (this->decoding_level_ > 2) {
          FastIDCT2D(this->current_image_->GetData(component_number - 1),
                     start_line + 8 * v_block, start_column + 8 * h_block,
                     line_length);
        }
      }
    }
    // If required, perform the color space tranform on the newly decoded
    // blocks.
  }
}

/**
 * \fn bool JPEGDecoder::IsMarker()
 *
 * \brief Tells whether the next two bytes in the stream are bytes of a marker.
 * This function does not increment the index count for the stream processed.
 *
 * \return True if the next bytes are a marker, else False.
 */
bool JPEGDecoder::IsMarker() {
  if (this->current_file_content_[this->current_index_] == 0xFF) {
    if (this->current_file_content_[this->current_index_ + 1] == 0x00) {
      return false;
    }
    return true;
  } else {
    return false;
  }
}

/**
 * \fn std::unique_ptr<unsigned char> JPEGDecoder::GetMarker()
 * \brief Returns the current marker if any, else throws a runtime error.
 *
 * All of the marker are of the form 0xFF.. with the .. being anything but 00.
 * The 00 are used to differanciate between the marker and values in the
 * encoded stream.
 *
 * This function expect to get a marker, it will throw an error if none found.
 * This function increments the index count for the stream processed.
 */
std::unique_ptr<unsigned char> JPEGDecoder::GetMarker() {
  std::stringstream error;
  unsigned char *marker;

  if (this->current_file_content_[this->current_index_] == 0xFF) {
    marker = new unsigned char[1];
    std::memcpy(marker,
                &(this->current_file_content_[this->current_index_ + 1]), 1);
    this->current_index_ = this->current_index_ + 2;
    return std::unique_ptr<unsigned char>(marker);
  } else {
    BOOST_LOG_TRIVIAL(error)
        << "Error while reading marker, 0xFF expected, but " << std::hex
        << (int)this->current_file_content_[this->current_index_]
        << " found at index: " << this->current_index_;
    throw std::runtime_error(error.str());
  }
}

/**
 * \fn void JPEGDecoder::Dequantize(const int component_number, const int
 * start_row, const int start_col, const QuantizationTable &table)
 *
 * \brief Dequantize a block of data.
 *
 * \param[in] component_number The component holding the data to dequantize.
 *
 * \param[in] start_row Cursor to the row at the start of the block.
 *
 * \param[in] start_col Cursor to the column at the start of the block.
 *
 * \param[in, out] table The quantization table to use.
 */
void JPEGDecoder::Dequantize(const int component_number, const int start_row,
                             const int start_col,
                             const QuantizationTable &table) {
  for (size_t row = 0; row < 8; row++) {
    for (size_t col = 0; col < 8; col++) {
#ifdef DEBUG
      this->current_image_->at(start_row + row, start_col + col,
                               component_number) *=
          table.qks_.at(ZZ_order[row * 8 + col]);
#else
      this->current_image_->at(start_row + row, start_col + col,
                               component_number) *=
          table.qks_[ZZ_order[row * 8 + col]];
#endif
    }
  }
}

/**
 * \fn void JPEGDecoder::InitializeLogger()
 *
 * \brief Function to initialize the logger.
 */
void JPEGDecoder::InitializeLogger() {
  switch (this->logging_level_) {
    case 0:
      boost::log::core::get()->set_filter(boost::log::trivial::severity >
                                          boost::log::trivial::fatal);
      break;
    case 1:
      boost::log::core::get()->set_filter(boost::log::trivial::severity ==
                                          boost::log::trivial::info);
      break;
    case 2:
      boost::log::core::get()->set_filter(boost::log::trivial::severity ==
                                          boost::log::trivial::debug);
      break;
    case 3:
      boost::log::core::get()->set_filter(boost::log::trivial::severity <=
                                          boost::log::trivial::info);
      break;
    case 4:
      boost::log::core::get()->set_filter(boost::log::trivial::severity <=
                                          boost::log::trivial::fatal);
      break;
    default:
      break;
  }
}

/**
 * \fn unsigned int JPEGDecoder::getImageSizeX()
 *
 * \brief Return the size of the image on the x axis (size of a row).
 *
 * \return The size of the image on the x axis (size of a row).
 */
unsigned int JPEGDecoder::getImageSizeX() {
  return this->frame_header_.number_of_samples_per_line_;
}

/**
 * \fn unsigned int JPEGDecoder::getImageSizeY()
 *
 * \brief
 */
unsigned int JPEGDecoder::getImageSizeY() {
  return this->frame_header_.number_of_lines_;
}

/**
 * \fn void JPEGDecoder::InitializeLogger()
 *
 * \brief Return the size of the image on the y axis (size of a col).
 *
 * \return The size of the image on the y axis (size of a col).
 */
int JPEGDecoder::getChannels() { return 3; }

/**
 * \fn int JPEGDecoder::getBlockPerLine()
 *
 * \brief Return the number of blocks per line.
 *
 * \return The number of blocks per line.
 */
int JPEGDecoder::getBlockPerLine() { return this->number_of_blocks_per_line; }

/**
 * \fn int JPEGDecoder::getBlockPerColumn()
 *
 * \brief Return the number of blocks per column.
 *
 * \return The number of blocks per column.
 */
int JPEGDecoder::getBlockPerColumn() {
  return this->number_of_blocks_per_column;
}
