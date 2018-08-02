#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"

/**
 * \fn cv::Mat JPEGDecoder::Decode(std::string filename, int level)
 * \brief Take a JPEG compress file as entry and output the decoded matrix.
 * 
 * \param[in] filename The file to be decoded.
 * \param[in] level The required level of decoding. 
 */
cv::Mat JPEGDecoder::Decode(std::string filename, int level)
{
  std::ifstream file_to_decode;
  int size, current_index = 0;
  unsigned char *file_content, *marker;

  // If the filename is the same, we assume to have the same image.
  if (filename.compare(this->current_filename_) == 0)
  {
    return this->current_image_.clone();
  }

  this->InitializeDecoder();
  this->decoding_level_ = level;

  file_to_decode.open(filename, std::ios::binary);

  if (file_to_decode.is_open())
  {
    // Getting the file content to be read.
    file_to_decode.seekg(0, std::ios::end);
    size = file_to_decode.tellg();
    file_to_decode.seekg(0, std::ios::beg);
    file_content = new unsigned char[size + 1];
    file_to_decode.read(reinterpret_cast<char *>(file_content), size);

    if (*(this->GetMarker(file_content, &current_index, 2)) == START_OF_IMAGE)
    {
      this->DecoderSetup();
      while (file_content[current_index] != END_OF_IMAGE)
      {
        marker = this->GetMarker(file_content, &current_index, 2);

        switch (*marker)
        {
        case APPO:
          this->ParseJFIFSegment(file_content, &current_index);
          break;
        case COMMENT:
          this->ParseComment(file_content, &current_index);
          break;
        case DEFINE_RESTART_INTERVAL:
          //TODO
          break;
        case DEFINE_QUANTIZATION_TABLE:
          this->ParseQuantizationTable(file_content, &current_index);
          std::cout << "quantization table parsed" << std::endl;
          break;
        case START_OF_FRAME_BASELINE:
          this->DecodeFrame(file_content, &current_index, FRAME_TYPE_BASELINE_DTC);
          std::cout << "frame header parsed" << std::endl;
          break;
        case START_OF_FRAME_PROGRESSIVE:
          this->DecodeFrame(file_content, &current_index, FRAME_TYPE_PROGRESSIVE);
          std::cout << "frame header parsed" << std::endl;
          break;
        case DEFINE_HUFFMAN_TABLE:
          this->ParseHuffmanTableSpecification(file_content, &current_index);
          std::cout << "huffman table specification parsed" << std::endl;
          break;
        default:
          std::cout << "I did not know how to parse the block : " << std::hex
                    << int(*marker) << std::endl;
          throw std::runtime_error("Error while processing the jpeg file.");
          break;
        }
      }
    }
    else
    {
      throw std::runtime_error("Expected SOI marker, but sosmething else found, cannot parse the file.");
    }
  }
  else
  {
    throw FileNotFoundException(filename);
  }

  return cv::Mat(300, 300, 25);
}

/**
 * \fn void JPEGDecoder::InitializeDecoder()
 * \brief Initialize all the values in the decoder for a new image to decode.
 */
void JPEGDecoder::InitializeDecoder()
{
  // We remove all elements from the maps of tables if any.
  this->quantization_tables_.clear();
  this->dc_huffman_tables_.clear();
  this->ac_huffman_tables_.clear();

  // We reset all the current header if required.
}

void JPEGDecoder::DecoderSetup()
{
  this->restart_interval = 0;
}

void JPEGDecoder::ParseJFIFSegment(unsigned char *file_content, int *index)
{
  unsigned char *marker;
  this->current_jfif_header = JFIFHeader();

  // Getting JFIF marker
  marker = this->GetMarker(file_content, index, 2, false);
  marker = this->GetMarker(file_content, index, 5, false);
  if (memcmp(marker, JFIF, 5) != 0)
  {
    throw std::runtime_error("Error while getting the JFIF marker.");
  }

  // Getting the version of the encoding.
  marker = this->GetMarker(file_content, index, 2, false);
  this->current_jfif_header.current_version_ = int(marker[0] << 8 | marker[1]);

  // Getting the units.
  marker = this->GetMarker(file_content, index, 1, false);
  this->current_jfif_header.current_unit_ = int(marker[0]);

  // Getting the Horizontal Pixel Density
  marker = this->GetMarker(file_content, index, 2, false);
  this->current_jfif_header.horizontal_pixel_density_ = int(marker[0] << 8 | marker[1]);

  // Getting the Vertical Pixel Density
  marker = this->GetMarker(file_content, index, 2, false);
  this->current_jfif_header.vertical_pixel_density_ = int(marker[0] << 8 | marker[1]);

  // Getting the Thumbnail Horizontal Pixel Count
  marker = this->GetMarker(file_content, index, 1, false);
  this->current_jfif_header.thumbnail_horizontal_pixel_count_ = int(marker[0]);

  // Getting the Thumbnail Vertical Pixel Count
  marker = this->GetMarker(file_content, index, 1, false);
  this->current_jfif_header.thumbnail_vertical_pixel_count_ = int(marker[0]);

  // If we have thumbnail count, we get the image.
  if (this->current_jfif_header.thumbnail_horizontal_pixel_count_ != 0 &&
      this->current_jfif_header.thumbnail_vertical_pixel_count_ != 0)
  {
    std::cout << "TODO" << std::endl;
  }
}

unsigned char *JPEGDecoder::GetMarker(unsigned char *file_content, int *index,
                                      int size, bool FF_expected)
{
  std::stringstream error;
  unsigned char *marker;

  if (FF_expected)
  {
    marker = new unsigned char[size - 1];
    std::memcpy(marker, &file_content[*index + 1], size - 1);

    if (file_content[*index] == 0xFF)
    {
      *index = *index + 2;
      return marker;
    }
    else
    {
      error << "Error while reading marker, 0xFF expected, but " << std::hex
            << (int)file_content[*index] << " found.";
      throw std::runtime_error(error.str());
    }
  }
  marker = new unsigned char[size];
  std::memcpy(marker, &file_content[*index], size);
  *index += size;
  return marker;
}

std::ostream &operator<<(std::ostream &out, const JPEGDecoder &decoder)
{
  out << "Format: JFIF\n";
  out << "Version : " << decoder.current_jfif_header.current_version_ << "\n";

  switch (decoder.current_jfif_header.current_unit_)
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

  out << "Horizontal pixel density : " << decoder.current_jfif_header.horizontal_pixel_density_
      << "\n";
  out << "Vertical pixel density : " << decoder.current_jfif_header.vertical_pixel_density_ << "\n";
  out << "Thumbnail horizontal pixel count : "
      << decoder.current_jfif_header.thumbnail_horizontal_pixel_count_ << "\n";
  out << "Thumbnail vertical pixel count : "
      << decoder.current_jfif_header.thumbnail_vertical_pixel_count_ << "\n";
  return out;
}

void JPEGDecoder::ParseQuantizationTable(unsigned char *file_content,
                                         int *index)
{
  unsigned int lq, temp;
  unsigned char pq, tq, mask_pq = 240, mask_tq = 15;
  unsigned int qk;
  QuantizationTable table_being_parsed = {}, insert_res;

  lq = int(file_content[*index] << 8 | file_content[*index + 1]);
  lq -= 2;
  *index += 2;
  while (lq > 0)
  {
    pq = file_content[*index];
    tq = (pq & mask_tq) > 4;
    pq = pq & mask_pq;

    table_being_parsed.pq_ = pq;

    *index += 1;
    if (pq == 1)
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[*index] << 8 | file_content[*index + 1]);
        *index += 2;
        table_being_parsed.qks_.push_back(qk);
      }
    }
    else
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[*index]);
        *index += 1;
        table_being_parsed.qks_.push_back(qk);
      }
    }
    if (!this->quantization_tables_
             .insert(std::pair<unsigned char, QuantizationTable>(
                 tq, table_being_parsed))
             .second)
    {
      this->quantization_tables_.erase(tq);
      this->quantization_tables_.insert(
          std::pair<unsigned char, QuantizationTable>(tq,
                                                      table_being_parsed));
    }
    temp = lq;
    lq -= (65 + 64 * pq);
    if (lq > temp)
    {
      lq = 0;
    }
  }
}

void JPEGDecoder::DecodeFrame(unsigned char *file_content, int *index, unsigned char encoding_process_type)
{
  unsigned char *marker;

  this->InterpretFrameHeader(file_content, index, encoding_process_type);

  do
  {
    marker = this->GetMarker(file_content, index, 2);
    if (*marker == START_OF_SCAN)
    {
      this->DecodeScan(file_content, index, encoding_process_type);
    }
    else
    {
      switch (*marker)
      {
      case COMMENT:
        this->ParseComment(file_content, index);
        break;
      case DEFINE_RESTART_INTERVAL:
        //TODO
        break;
      case DEFINE_QUANTIZATION_TABLE:
        this->ParseQuantizationTable(file_content, index);
        std::cout << "quantization table parsed" << std::endl;
        break;
      case DEFINE_HUFFMAN_TABLE:
        this->ParseHuffmanTableSpecification(file_content, index);
        std::cout << "huffman table specification parsed" << std::endl;
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

  // Now that the whole frame is huffman decoded, depending on the selected level, we keep extracting or stop here

  switch (this->decoding_level_)
  {
  case 1:
    // One is huffman level, not de-quantized, so we do nothing.
    break;
  case 2:
    // Two is de-quantized huffman decoded values.
    break;
  case 3:
    // Three is IDCT coded image.
    break;
  case 4:
    // Four is YCbCr image.
    break;
  case 5:
    // Five is normal image.
  default:
    throw std::runtime_error("Unexpected level of decoding, exiting.");
    break;
  }
}

void JPEGDecoder::DecodeScan(unsigned char *file_content, int *index, unsigned char encoding_process_type)
{
  this->ParseScanHeader(file_content, index);
  unsigned int m = 0;

  do
  {
    m += 1;

    switch (encoding_process_type)
    {
    case FRAME_TYPE_BASELINE_DTC:
      this->DecodeRestartIntervalBaseline(file_content, index);
      break;

    default:
      break;
    }
  } while (m < this->restart_interval);
}

void JPEGDecoder::DecodeRestartIntervalBaseline(unsigned char *file_content, int *index)
{
  int n = 0;
  unsigned char decoded_dc, diff;
  this->ResetDecoderBaseline();
  int component_number = 0;
  unsigned char dc_table_index, ac_table_index;
  cv::Mat new_block = cv::Mat(8, 8, 1);

  while (!this->IsMarker(file_content, *index))
  {
    while (n < this->data_unit_per_mcu_)
    {
      if (this->current_frame_header_.number_image_component > 1)
      {
        // Then we are interleaved.
        dc_table_index = this->current_scan_.scan_components_specification_parameters_.at(component_number).first;
        ac_table_index = this->current_scan_.scan_components_specification_parameters_.at(component_number).second;

        //We decode the DC component.
        decoded_dc = this->DecodeBaseline(file_content, index, this->dc_huffman_tables_.at(dc_table_index));
        diff = this->ReceiveBaseline(decoded_dc);
        diff = this->ExtendedBaseline(diff, decoded_dc);

        new_block.at<uchar>(0, 0, 0) = diff;

        // We decode the ac components.
        this->DecodeACCoefficients(file_content, index, &new_block, this->ac_huffman_tables_.at(ac_table_index));

        if (component_number == this->current_frame_header_.number_image_component - 1)
        {
          component_number = 0;
        }
        else
        {
          component_number += 1;
        }
      }
    }
  }
}

int JPEGDecoder::ReceiveBaseline(unsigned char decoded_dc)
{
  int value = 0, i = 0;

  while (i != decoded_dc)
  {
    i += i;
    value = (value << 1) + this->NextBit(this->current_file_content_, &i);
  }

  return value;
}

int JPEGDecoder::ExtendedBaseline(unsigned char diff, unsigned char decoded_dc)
{
  int value = 1;

  value = value << (decoded_dc - 1);
  if (diff < value)
  {
    value = (-1 << decoded_dc) + 1;
    diff += value;
  }
  return value;
}

void JPEGDecoder::DecodeRestartIntervalProgressive(unsigned char *file_content, int *index)
{
  bool mcu = true; // To define
  this->ResetDecoderProgressive();

  while (!this->IsMarker(file_content, *index))
  {
    this->DecodeMCUProgressive(file_content, index);
  }
}

void JPEGDecoder::DecodeMCUProgressive(unsigned char *file_content, int *index)
{
  int n = 0;

  while (n < this->data_unit_per_mcu_)
  {
    // ADD the decoding function.
  }
}

void JPEGDecoder::InterpretFrameHeader(unsigned char *file_content, int *index, unsigned char encoding_process_type)
{
  unsigned int header_length;
  unsigned char number_of_image_component_in_frame, mask_h = 240, mask_v = 15, c, horizontal_sampling_factor, vertical_sampling_factor, quantization_table_destination_selector;

  this->current_frame_header_.encoding_process_type_ = encoding_process_type;

  header_length = int(file_content[*index] << 8 | file_content[*index + 1]);
  this->current_frame_header_.sample_precision_ = file_content[*index + 2];
  this->current_frame_header_.number_of_lines_ = int(file_content[*index + 3] << 8 | file_content[*index + 4]);
  this->current_frame_header_.number_of_samples_per_line_ = int(file_content[*index + 5] << 8 | file_content[*index + 6]);
  number_of_image_component_in_frame = file_content[*index + 7];
  *index += 8;

  for (size_t i = 0; i < number_of_image_component_in_frame; i++)
  {
    std::vector<unsigned char> components_vector;
    c = file_content[*index];
    horizontal_sampling_factor = (file_content[*index + 1] & mask_h) >> 4;
    vertical_sampling_factor = file_content[*index + 1] & mask_v;
    quantization_table_destination_selector = file_content[*index + 2];
    components_vector.push_back(horizontal_sampling_factor);
    components_vector.push_back(vertical_sampling_factor);
    components_vector.push_back(quantization_table_destination_selector);
    this->current_frame_header_.component_signification_parameters_.insert(std::pair<unsigned char, std::vector<unsigned char>>(c, components_vector));
    *index += 3;
  }
}

void JPEGDecoder::ParseHuffmanTableSpecification(unsigned char *file_content, int *index)
{
  unsigned int table_definition_length, counter, temp, sum_code_length = 0;
  unsigned char table_destination_identifier, table_class, mask_high = 240, mask_low = 15, number_of_huffman_codes;

  HuffmanTable table_being_parsed = {};

  table_definition_length = int(file_content[*index] << 8 | file_content[*index + 1]);
  table_definition_length -= 2;
  *index += 2;

  while (table_definition_length > 0)
  {
    sum_code_length = 0;
    table_class = (file_content[*index] & mask_high) >> 4;
    table_destination_identifier = file_content[*index] & mask_low;

    *index += 1;

    for (size_t i = 0; i < 16; i++)
    {
      number_of_huffman_codes = file_content[*index];
      table_being_parsed.bits.push_back(number_of_huffman_codes);
      sum_code_length += number_of_huffman_codes;
      *index += 1;
    }

    for (size_t i = 0; i < 16; i++)
    {
      counter = table_being_parsed.bits.at(i);
      table_being_parsed.huffvals.push_back(std::vector<unsigned char>());

      for (size_t j = 0; j < counter; j++)
      {
        table_being_parsed.huffvals.at(i).push_back(file_content[*index]);
        *index += 1;
      }
    }
    temp = table_definition_length;
    table_definition_length = table_definition_length - 17 - sum_code_length;
    if (temp < table_definition_length)
    {
      table_definition_length = 0;
    }
  }

  // Parsing the retrieved values.
  table_being_parsed.huffsize = this->GenerateSizeTable(table_being_parsed.bits);
  this->GenerateCodeTable(&table_being_parsed);

  // If 0, DC table, else AC.
  if (table_class == 0)
  {
    if (!this->dc_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed)).second)
    {
      this->dc_huffman_tables_.erase(table_destination_identifier);
      this->dc_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed));
    }
  }
  else
  {
    if (!this->ac_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed)).second)
    {
      this->ac_huffman_tables_.erase(table_destination_identifier);
      this->ac_huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed));
    }
  }
}

void JPEGDecoder::ParseScanHeader(unsigned char *file_content, int *index)
{
  // unsigned int scan_header_length;
  unsigned char image_component_in_scan, scan_component_selector, dc_table_destination_selector, ac_table_destination_selector, mask_high = 240, mask_low = 15;
  ScanHeader current_scan = {};

  // To do Check on the length of the header.
  // scan_header_length = int(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;
  image_component_in_scan = file_content[*index];
  *index += 1;

  this->current_scan_.number_of_image_components_ = image_component_in_scan;

  for (size_t i = 0; i < image_component_in_scan; i++)
  {
    scan_component_selector = file_content[*index];
    dc_table_destination_selector = (file_content[*index + 1] & mask_high) >> 4;
    ac_table_destination_selector = file_content[*index + 1] & mask_low;
    *index += 2;
    if (!this->current_scan_.scan_components_specification_parameters_
             .insert(std::pair<unsigned char, std::pair<unsigned char, unsigned char>>(
                 scan_component_selector, std::pair<unsigned char, unsigned char>(dc_table_destination_selector, ac_table_destination_selector)))
             .second)
    {
      this->quantization_tables_.erase(scan_component_selector);
      this->current_scan_.scan_components_specification_parameters_
          .insert(std::pair<unsigned char, std::pair<unsigned char, unsigned char>>(
              scan_component_selector, std::pair<unsigned char, unsigned char>(dc_table_destination_selector, ac_table_destination_selector)));
    }
  }

  this->current_scan_.start_of_spectral_selection_ = file_content[*index];
  this->current_scan_.end_of_spectral_selection_ = file_content[*index + 1];
  this->current_scan_.approximation_high_bit_ = (file_content[*index + 2] & mask_high) >> 4;
  this->current_scan_.approximation_low_bit_ = file_content[*index + 2] & mask_low;
}

void JPEGDecoder::ParseComment(unsigned char *file_content, int *index)
{
  unsigned int comment_length;

  comment_length = (unsigned int)(file_content[*index] << 8 | file_content[*index + 1]);

  std::cout << "Comment in the file : ";
  for (size_t i = 0; i < comment_length - 2; i++)
  {
    std::cout << file_content[*index + 2 + i];
  }

  *index += comment_length;
}

unsigned char JPEGDecoder::NextBit(unsigned char *file_content, int *index)
{
  char current_byte = file_content[*index], bit;

  if (this->next_bit_count_ == 0)
  {
    this->next_bit_count_ = 8;
    if (current_byte == 0xFF)
    {
      if (!(file_content[*index + 1] == 0x00))
      {
        if (file_content[*index + 1] == DEFINE_NUMBER_OF_LINE)
        {
          // TODO //
          // Should process the dnl marker
        }
        else
        {
          throw std::runtime_error("Error while parsing the file, DNL byte expected but something else found");
        }
      }
      else
      {
        if ((this->next_bit_count_ - 1) == 0)
        {
          *index += 1;
        }
      }
    }
  }

  bit = current_byte << (8 - this->next_bit_count_);
  bit = bit >> 7;
  this->next_bit_count_ -= 1;

  if (this->next_bit_count_ == 0)
  {
    *index += 1;
  }
  return bit;
}

unsigned char JPEGDecoder::DecodeBaseline(unsigned char *file_content, int *index, HuffmanTable used_table)
{
  int i = 1, j;
  char code;

  code = this->NextBit(file_content, index);
  while (code > used_table.max_code.at(i))
  {
    i += 1;
    code = (code << 1) + this->NextBit(file_content, index);
  }

  j = used_table.val_pointer.at(j);
  j = j + code - used_table.max_code.at(i);
  return used_table.huffvals.at(i).at(j);
}

std::vector<unsigned char> JPEGDecoder::GenerateSizeTable(std::vector<unsigned char> bits)
{
  unsigned char k = 0, i = 1, j = 1;
  std::vector<unsigned char> huffsize;

  while (i < 16)
  {
    if (j > bits.at(i))
    {
      i += 1;
      j = 1;
    }
    else
    {
      huffsize.push_back(i);
      k += 1;
      j += 1;
    }
  }
  huffsize.push_back(i);
  this->last_k_ = k;
}

void JPEGDecoder::GenerateCodeTable(HuffmanTable *table_to_fill)
{
  bool smallest_processed;
  int k = 0, code = 0;
  unsigned char si = table_to_fill->huffsize.at(0);
  std::vector<int> huffcode(table_to_fill->huffsize.size(), 0), min_code(16, 0), max_code(16, 0), val_pointer(16, 0);

  while (true)
  {
    smallest_processed = false;
    do
    {
      if (!smallest_processed)
      {
        val_pointer.at(si) = k;
        min_code.at(si) = code;
        smallest_processed = true;
      }
      huffcode.at(k) = code;
      code += 1;
      k += 1;
    } while (table_to_fill->huffsize.at(k) == si);
    max_code.at(si) = code;
    if (table_to_fill->huffsize.at(k) == 0)
    {
      table_to_fill->huffcode = huffcode;
      table_to_fill->max_code = max_code;
      table_to_fill->min_code = min_code;
      table_to_fill->val_pointer = val_pointer;
      return;
    }

    do
    {
      code = code << 1;
      si = si + 1;
    } while (table_to_fill->huffsize.at(k) != si);
  }
}

void JPEGDecoder::ResetDecoderBaseline()
{
}

bool JPEGDecoder::IsMarker(unsigned char *file_content, int index)
{
  return true;
}

/**
 * \fn void JPEGDecoder::DecodeACCoefficients(unsigned char *file_content, int *index, cv::Mat *new_block)
 * \brief Decode the AC coefficients for the baseline procedure and store them in the new_block.
 * 
 * \param[in] file_content Pointer to the file being parsed.
 * \param[in, out] index, Pointer to the position of the cursor in the file being parsed.
 * \param[in, out] new_block Pointer to the block of data to contain the decoded values.
 */
void JPEGDecoder::DecodeACCoefficients(unsigned char *file_content, int *index, cv::Mat *new_block, HuffmanTable used_table)
{
  unsigned char k = 1, rs = 0, ssss = 0, rrrr = 0, r = 0;
  std::vector<unsigned char> ZZ(63, 0);
  bool out_condition = false;

  do
  {
    rs = this->DecodeBaseline(file_content, index, used_table);
    ssss = rs % 16;
    rrrr = rs >> rs;
    r = rrrr;
    if (ssss == 0)
    {
      if (r != 15)
      {
        k = k + 16;
      }
      else
      {
        out_condition = true;
      }
    }
    else
    {
      k = k + r;
      ZZ.at(k) = this->DecodeZZ(k, ssss);
      if (k == 63)
      {
        out_condition = true;
      }
    }
  } while (!out_condition);
}

void JPEGDecoder::ResetDecoderProgressive() {}

/**
 * \fn unsigned char DecodeZZ(unsigned char k, unsigned char ssss)
 * \brief Decode the coefficient in the zigzag order.
 * 
 * \param[in] k The position of the coefficient.
 * \param[in] ssss the low bit of ?
 */
unsigned char JPEGDecoder::DecodeZZ(unsigned char k, unsigned char ssss)
{
  unsigned char return_value;
  return_value = this->ReceiveBaseline(ssss);
  return_value = this->ExtendedBaseline(return_value, ssss);
}