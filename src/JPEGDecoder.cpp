#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "JPEGDecoder.hpp"
#include "JPEGException.hpp"

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

  file_to_decode.open(filename, std::ios::binary);

  if (file_to_decode.is_open())
  {
    file_to_decode.seekg(0, std::ios::end);
    size = file_to_decode.tellg();
    file_to_decode.seekg(0, std::ios::beg);
    file_content = new unsigned char[size + 1];
    file_to_decode.read(reinterpret_cast<char *>(file_content), size);

    if (this->GetFileInformation(file_content, &current_index))
    {
      while (file_content[current_index] != END_OF_IMAGE)
      {
        std::cout << current_index << std::endl;
        marker = this->GetMarker(file_content, &current_index, 2);

        switch (*marker)
        {
        case COMMENT:
          this->ParseComment(file_content, &current_index);
          break;
        case DEFINE_QUANTIZATION_TABLE:
          this->ParseQuantizationTable(file_content, &current_index);
          std::cout << "quantization table parsed" << std::endl;
          break;
        case START_OF_FRAME_BASELINE:
          this->ParseFrameHeader(file_content, &current_index, FRAME_TYPE_BASELINE_DTC);
          std::cout << "frame header parsed" << std::endl;
          break;
        case START_OF_FRAME_PROGRESSIVE:
          this->ParseFrameHeader(file_content, &current_index, FRAME_TYPE_PROGRESSIVE);
          std::cout << "frame header parsed" << std::endl;
          break;
        case DEFINE_HUFFMAN_TABLE:
          this->ParseHuffmanTableSpecification(file_content, &current_index);
          std::cout << "huffman table specification parsed" << std::endl;
          break;
        case START_OF_SCAN:
          this->ParseScanHeader(file_content, &current_index);
          std::cout << "Scan header parsed" << std::endl;
          break;
        case END_OF_IMAGE:
          std::cout << "end of image." << std::endl;
          break;
        default:
          std::cout << "I did not know how to parse the block : " << std::hex
                    << int(*marker) << std::endl;
          throw std::runtime_error("Error while processing the jpeg file.");
          break;
        }
      }

      std::cout << "Done parsing the image file, we will start decoding the "
                   "image according to the level specified by the user."
                << std::endl;

      switch (level)
      {
      case 1:
        std::cout << "We output the image as a encoded before the DCT, but "
                     "after over-sampling."
                  << std::endl;
        break;
      case 2:
        std::cout << "We output the image as a encoded after the DCT."
                  << std::endl;
        break;
      case 3:
        std::cout << "We output the image as encoded before the huffman "
                     "algorithm, but after quantization."
                  << std::endl;
        break;
      case 4:
        std::cout
            << "We output the image as encoded with the huffman algorithm."
            << std::endl;
        break;
      case 5:
        std::cout
            << "We output the image as a encoded with the huffman algorithm."
            << std::endl;
        break;
      default:
        break;
      }
    }
    else
    {
      throw std::runtime_error(
          "Error while getting the information on the file to parse.");
    }
  }
  else
  {
    throw FileNotFoundException(filename);
  }

  return cv::Mat(300, 300, 25);
}

bool JPEGDecoder::GetFileInformation(unsigned char *file_content, int *index)
{
  unsigned char *marker;

  // Start of image marker.
  marker = this->GetMarker(file_content, index, 2);
  if (*marker != START_OF_IMAGE)
  {
    return false;
  }

  // Getting APP0 marker.
  marker = this->GetMarker(file_content, index, 2);

  if (*marker != APPO)
  {
    return false;
  }

  // Getting JFIF marker
  marker = this->GetMarker(file_content, index, 2, false);
  marker = this->GetMarker(file_content, index, 5, false);
  if (memcmp(marker, JFIF, 5) != 0)
  {
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
  if (this->thumbnail_horizontal_pixel_count_ != 0 &&
      this->thumbnail_vertical_pixel_count_ != 0)
  {
    std::cout << "TODO" << std::endl;
  }

  return true;
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

  out << "Horizontal pixel density : " << decoder.horizontal_pixel_density_
      << "\n";
  out << "Vertical pixel density : " << decoder.vertical_pixel_density_ << "\n";
  out << "Thumbnail horizontal pixel count : "
      << decoder.thumbnail_horizontal_pixel_count_ << "\n";
  out << "Thumbnail vertical pixel count : "
      << decoder.thumbnail_vertical_pixel_count_ << "\n";
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

void JPEGDecoder::ParseFrameHeader(unsigned char *file_content, int *index, unsigned char encoding_process_type)
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
  unsigned char table_destination_identifier, mask_high = 240, mask_low = 15, number_of_huffman_codes;

  HuffmanTable table_being_parsed = {};

  table_definition_length = int(file_content[*index] << 8 | file_content[*index + 1]);
  table_definition_length -= 2;
  *index += 2;

  while (table_definition_length > 0)
  {
    sum_code_length = 0;
    table_being_parsed.table_class_ = (file_content[*index] & mask_high) >> 4;
    table_destination_identifier = file_content[*index] & mask_low;

    *index += 1;

    for (size_t i = 0; i < 16; i++)
    {
      number_of_huffman_codes = file_content[*index];
      table_being_parsed.number_of_codes_of_length_i_.push_back(number_of_huffman_codes);
      sum_code_length += number_of_huffman_codes;
      *index += 1;
    }

    for (size_t i = 0; i < 16; i++)
    {
      counter = table_being_parsed.number_of_codes_of_length_i_.at(i);
      table_being_parsed.huffman_code_associated_values_.push_back(std::vector<unsigned int>());

      for (size_t j = 0; j < counter; j++)
      {
        table_being_parsed.huffman_code_associated_values_.at(i).push_back(file_content[*index]);
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
  if (!this->huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed)).second)
  {
    this->quantization_tables_.erase(table_destination_identifier);
    this->huffman_tables_.insert(std::pair<unsigned char, HuffmanTable>(table_destination_identifier, table_being_parsed));
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

char JPEGDecoder::NextBit(unsigned char *file_content, int *index)
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

char JPEGDecoder::Decode(unsigned char *file_content, int *index, HuffmanTable used_table)
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
  return used_table.huffman_code_associated_values_.at(j);
}

std::vector<char> JPEGDecoder::GenerateSizeTable(std::vector<char> bits)
{
  int k = 0, i = 1, j = 1;
  std::vector<char> huffsize;

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

std::vector<int> JPEGDecoder::GenerateCodeTable(std::vector<char> huffsize)
{
  bool smallest_processed;
  int k = 0, code = 0;
  char si = huffsize.at(0);
  std::vector<int> huffcode(huffsize.size(), 0), min_code(16, 0), max_code(16, 0), val_pointer(16, 0);

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
    } while (huffsize.at(k) == si);
    max_code.at(si) = code;
    if (huffsize.at(k) == 0)
    {
      return huffcode;
    }

    do
    {
      code = code << 1;
      si = si + 1;
    } while (huffsize.at(k) != si);
  }
}