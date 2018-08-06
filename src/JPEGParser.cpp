#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "JPEGParser.hpp"
#include "JPEGType.hpp"

/**
 * \fn JFIFHeader ParseJFIFSegment(unsigned char* file_content, int*
 * index) \brief Parse the JFIF Segment at the begining of most of the jpeg
 * files.
 */
JFIFHeader ParseJFIFSegment(unsigned char* file_content, unsigned int* index) {
  unsigned char* marker;
  JFIFHeader jfif_header = JFIFHeader();
  unsigned int length;

  // Getting the length of the segment
  length = (unsigned int)(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;
  if (memcmp(&(file_content[*index]), JFIF, 5) != 0) {
    throw std::runtime_error("Error while getting the JFIF marker.");
  }
  *index += 5;

  // Getting the version of the encoding.
  jfif_header.current_version_ =
      int(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;

  // Getting the units.
  jfif_header.current_unit_ = int(file_content[*index]);
  *index += 1;

  // Getting the Horizontal Pixel Density
  jfif_header.horizontal_pixel_density_ =
      int(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;

  // Getting the Vertical Pixel Density
  jfif_header.vertical_pixel_density_ =
      int(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;

  // Getting the Thumbnail Horizontal Pixel Count
  jfif_header.thumbnail_horizontal_pixel_count_ = int(file_content[*index]);
  *index += 1;

  // Getting the Thumbnail Vertical Pixel Count
  jfif_header.thumbnail_vertical_pixel_count_ = int(file_content[*index]);
  *index += 1;

  // If we have thumbnail count, we get the image.
  if (jfif_header.thumbnail_horizontal_pixel_count_ != 0 &&
      jfif_header.thumbnail_vertical_pixel_count_ != 0) {
    std::cout << "TODO" << std::endl;
  }
  return jfif_header;
}

/**
 * std::string ParseComment(unsigned char* file_content, int* index)
 * Parse the comment in the comment block. A comment can be anything and is not
 * processed by the decoder.
 */
std::string ParseComment(unsigned char* file_content, unsigned int* index) {
  unsigned int comment_length;
  std::string comment;

  comment_length =
      (unsigned int)(file_content[*index] << 8 | file_content[*index + 1]);

  comment = std::string(file_content[*index + 2], comment_length - 2);

  *index += comment_length;
  return comment;
}

/**
 * \fn void ParseFrameHeader(unsigned char encoding_process_type)
 * \brief Parse the Frame header and store all the information in it.
 *
 * \param[in] encoding_process_type The type of encoding used.
 */
FrameHeader ParseFrameHeader(unsigned char* file_content, unsigned int* index,
                             unsigned char encoding_process_type) {
  FrameHeader frame_header = FrameHeader();

  // unsigned int header_length;
  unsigned char number_of_image_component_in_frame,
      mask_h = 240, mask_v = 15, c, horizontal_sampling_factor,
      vertical_sampling_factor, quantization_table_destination_selector;

  frame_header.encoding_process_type_ = encoding_process_type;

  // To use to check header size.
  // header_length = int(file_content[*index] << 8 |
  // file_content[*index + 1]);

  frame_header.sample_precision_ = file_content[*index + 2];

  frame_header.number_of_lines_ =
      int(file_content[*index + 3] << 8 | file_content[*index + 4]);

  frame_header.number_of_samples_per_line_ =
      int(file_content[*index + 5] << 8 | file_content[*index + 6]);

  frame_header.number_of_image_component = file_content[*index + 7];

  *index += 8;

  for (size_t i = 0; i < frame_header.number_of_image_component; i++) {
    std::vector<unsigned char> components_vector;

    c = file_content[*index];

    horizontal_sampling_factor = (file_content[*index + 1] & mask_h) >> 4;
    vertical_sampling_factor = file_content[*index + 1] & mask_v;

    quantization_table_destination_selector = file_content[*index + 2];

    components_vector.push_back(horizontal_sampling_factor);
    components_vector.push_back(vertical_sampling_factor);
    components_vector.push_back(quantization_table_destination_selector);

    frame_header.component_signification_parameters_.insert(
        std::pair<unsigned char, std::vector<unsigned char>>(
            c, components_vector));
    *index += 3;
  }
}

/**
 * \fn ScanHeader ParseScanHeader(unsigned char* file_content, int*
 * index) \brief Parse the scan header.
 */
ScanHeader ParseScanHeader(unsigned char* file_content, unsigned int* index) {
  ScanHeader scan = ScanHeader();

  // unsigned int scan_header_length;
  unsigned char image_component_in_scan, scan_component_selector,
      dc_table_destination_selector, ac_table_destination_selector,
      mask_high = 240, mask_low = 15;
  ScanHeader current_scan = {};

  // To do Check on the length of the header.
  // scan_header_length = int(file_content[*index] << 8 |
  // file_content[*index + 1]);
  *index += 2;

  image_component_in_scan = file_content[*index];
  *index += 1;

  scan.number_of_image_components_ = image_component_in_scan;

  for (size_t i = 0; i < image_component_in_scan; i++) {
    scan_component_selector = file_content[*index];
    dc_table_destination_selector = (file_content[*index + 1] & mask_high) >> 4;
    ac_table_destination_selector = file_content[*index + 1] & mask_low;
    *index += 2;
    if (!scan.scan_components_specification_parameters_
             .insert(std::pair<unsigned char,
                               std::pair<unsigned char, unsigned char>>(
                 scan_component_selector,
                 std::pair<unsigned char, unsigned char>(
                     dc_table_destination_selector,
                     ac_table_destination_selector)))
             .second) {
      scan.scan_components_specification_parameters_.insert(
          std::pair<unsigned char, std::pair<unsigned char, unsigned char>>(
              scan_component_selector, std::pair<unsigned char, unsigned char>(
                                           dc_table_destination_selector,
                                           ac_table_destination_selector)));
    }
  }

  scan.start_of_spectral_selection_ = file_content[*index];
  scan.end_of_spectral_selection_ = file_content[*index + 1];
  scan.approximation_high_bit_ = (file_content[*index + 2] & mask_high) >> 4;
  scan.approximation_low_bit_ = file_content[*index + 2] & mask_low;

  return scan;
}

/**
 * \fn QuantizationTable ParseQuantizationTable(unsigned char*
 * file_content, int* index)
 * \brief Parse the quantization table and store the
 * information at the specified location.
 */
QuantizationTable ParseQuantizationTable(unsigned char* file_content,
                                         unsigned int* index) {
  unsigned int lq, temp;
  unsigned char pq, tq, mask_pq = 240, mask_tq = 15;
  unsigned int qk;
  QuantizationTable table_being_parsed = {};

  lq = int(file_content[*index] << 8 | file_content[*index + 1]);
  lq -= 2;
  *index += 2;
  while (lq > 0) {
    pq = file_content[*index];
    tq = (pq & mask_tq) > 4;
    pq = pq & mask_pq;

    table_being_parsed.pq_ = pq;

    *index += 1;
    if (pq == 1) {
      for (size_t i = 0; i < 64; i++) {
        qk = int(file_content[*index] << 8 | file_content[*index + 1]);
        *index += 2;
        table_being_parsed.qks_.push_back(qk);
      }
    } else {
      for (size_t i = 0; i < 64; i++) {
        qk = int(file_content[*index]);
        *index += 1;
        table_being_parsed.qks_.push_back(qk);
      }
    }
    temp = lq;
    lq -= (65 + 64 * pq);
    if (lq > temp) {
      lq = 0;
    }
  }
  return table_being_parsed;
}

/**
 * \fn std::vector<std::pair<unsigned char, HuffmanTable>>
 * ParseHuffmanTableSpecification(unsigned char* file_content, int*
 * index) \brief Parse the huffman table and create the associated tables.
 */
std::vector<std::pair<unsigned char, HuffmanTable>>
ParseHuffmanTableSpecification(unsigned char* file_content,
                               unsigned int* index) {
  unsigned int table_definition_length, counter, temp, sum_code_length = 0;
  unsigned char table_destination_identifier, table_class,
      mask_high = 240, mask_low = 15, number_of_huffman_codes;
  std::vector<std::pair<unsigned char, HuffmanTable>> output;

  HuffmanTable table_being_parsed;

  table_definition_length =
      int(file_content[*index] << 8 | file_content[*index + 1]);
  table_definition_length -= 2;
  *index += 2;

  while (table_definition_length > 0) {
    sum_code_length = 0;
    table_class = (file_content[*index] & mask_high) >> 4;
    table_destination_identifier = file_content[*index] & mask_low;

    *index += 1;

    for (size_t i = 0; i < 16; i++) {
      number_of_huffman_codes = file_content[*index];
      table_being_parsed.bits.push_back(number_of_huffman_codes);
      sum_code_length += number_of_huffman_codes;
      *index += 1;
    }

    for (size_t i = 0; i < 16; i++) {
      counter = table_being_parsed.bits.at(i);

      for (size_t j = 0; j < counter; j++) {
        table_being_parsed.huffvals.push_back(file_content[*index]);
        *index += 1;
      }
    }
    temp = table_definition_length;
    table_definition_length = table_definition_length - 17 - sum_code_length;
    if (temp < table_definition_length) {
      table_definition_length = 0;
    }
  }

  /*   // Parsing the retrieved values.
    table_being_parsed.huffsize =
        this->huffman_decoder_.GenerateSizeTable(table_being_parsed.bits);
    table_being_parsed.huffcode =
        this->huffman_decoder_.GenerateCodeTable(table_being_parsed.huffsize);
    this->huffman_decoder_.DecoderTables(&table_being_parsed);

    // If 0, DC table, else AC.
    if (table_class == 0) {
      if (!this->dc_huffman_tables_
               .insert(std::pair<unsigned char, HuffmanTable>(
                   table_destination_identifier, table_being_parsed))
               .second) {
        this->dc_huffman_tables_.erase(table_destination_identifier);
        this->dc_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(
            table_destination_identifier, table_being_parsed));
      }
    } else {
      if (!this->ac_huffman_tables_
               .insert(std::pair<unsigned char, HuffmanTable>(
                   table_destination_identifier, table_being_parsed))
               .second) {
        this->ac_huffman_tables_.erase(table_destination_identifier);
        this->ac_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(
            table_destination_identifier, table_being_parsed));
      }
    } */
}
