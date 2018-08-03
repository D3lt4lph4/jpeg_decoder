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
  unsigned char *marker;

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
    this->current_file_content_ = new unsigned char[size + 1];
    file_to_decode.read(reinterpret_cast<char *>(this->current_file_content_), size);

    if (*(this->GetMarker()) == START_OF_IMAGE)
    {
      this->DecoderSetup();
      while (this->current_file_content_[current_index] != END_OF_IMAGE)
      {
        marker = this->GetMarker();

        switch (*marker)
        {
        case APPO:
          this->ParseJFIFSegment();
          break;
        case COMMENT:
          this->ParseComment();
          break;
        case DEFINE_RESTART_INTERVAL:
          //TODO
          break;
        case DEFINE_QUANTIZATION_TABLE:
          this->ParseQuantizationTable();
          break;
        case START_OF_FRAME_BASELINE:
          this->DecodeFrame(FRAME_TYPE_BASELINE_DTC);
          break;
        case START_OF_FRAME_PROGRESSIVE:
          this->DecodeFrame(FRAME_TYPE_PROGRESSIVE);
          break;
        case DEFINE_HUFFMAN_TABLE:
          this->ParseHuffmanTableSpecification();
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

void JPEGDecoder::ParseJFIFSegment()
{
  unsigned char *marker;
  this->current_jfif_header = JFIFHeader();

  // Getting JFIF marker
  marker = this->GetMarker();
  marker = this->GetMarker();
  if (memcmp(marker, JFIF, 5) != 0)
  {
    throw std::runtime_error("Error while getting the JFIF marker.");
  }

  // Getting the version of the encoding.
  marker = this->GetMarker();
  this->current_jfif_header.current_version_ = int(marker[0] << 8 | marker[1]);

  // Getting the units.
  marker = this->GetMarker();
  this->current_jfif_header.current_unit_ = int(marker[0]);

  // Getting the Horizontal Pixel Density
  marker = this->GetMarker();
  this->current_jfif_header.horizontal_pixel_density_ = int(marker[0] << 8 | marker[1]);

  // Getting the Vertical Pixel Density
  marker = this->GetMarker();
  this->current_jfif_header.vertical_pixel_density_ = int(marker[0] << 8 | marker[1]);

  // Getting the Thumbnail Horizontal Pixel Count
  marker = this->GetMarker();
  this->current_jfif_header.thumbnail_horizontal_pixel_count_ = int(marker[0]);

  // Getting the Thumbnail Vertical Pixel Count
  marker = this->GetMarker();
  this->current_jfif_header.thumbnail_vertical_pixel_count_ = int(marker[0]);

  // If we have thumbnail count, we get the image.
  if (this->current_jfif_header.thumbnail_horizontal_pixel_count_ != 0 &&
      this->current_jfif_header.thumbnail_vertical_pixel_count_ != 0)
  {
    std::cout << "TODO" << std::endl;
  }
}

void JPEGDecoder::ParseQuantizationTable()
{
  unsigned int lq, temp;
  unsigned char pq, tq, mask_pq = 240, mask_tq = 15;
  unsigned int qk;
  QuantizationTable table_being_parsed = {}, insert_res;

  lq = int(this->current_file_content_[this->current_index_] << 8 | this->current_file_content_[this->current_index_ + 1]);
  lq -= 2;
  this->current_index_ += 2;
  while (lq > 0)
  {
    pq = this->current_file_content_[this->current_index_];
    tq = (pq & mask_tq) > 4;
    pq = pq & mask_pq;

    table_being_parsed.pq_ = pq;

    this->current_index_ += 1;
    if (pq == 1)
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[this->current_index_] << 8 | file_content[this->current_index_ + 1]);
        this->current_index_ += 2;
        table_being_parsed.qks_.push_back(qk);
      }
    }
    else
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[this->current_index_]);
        this->current_index_ += 1;
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

void JPEGDecoder::DecodeFrame(unsigned char this->current_file_content_, int this->current_index_, unsigned char encoding_process_type)
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
        std::cout << "huffman table specification parsed_1" << std::endl;
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

void JPEGDecoder::DecodeScan(unsigned char this->current_file_content_, int this->current_index_, unsigned char encoding_process_type)
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

void JPEGDecoder::DecodeRestartIntervalBaseline(unsigned char this->current_file_content_, int this->current_index_)
{
  int n = 0;
  unsigned char decoded_dc, diff;
  this->ResetDecoderBaseline();
  unsigned int component_number = 1;
  unsigned char dc_table_index, ac_table_index;
  this->data_unit_per_mcu_ = 3;
  cv::Mat new_block = cv::Mat(8, 8, 1);

  while (!this->IsMarker(file_content, this->current_index_))
  {
    n = 0;
    while (n < this->data_unit_per_mcu_)
    {
      if (this->current_frame_header_.number_of_image_component > 1)
      {
        // Then we are interleaved.
        dc_table_index = this->current_scan_.scan_components_specification_parameters_.at((unsigned char)component_number).first;
        ac_table_index = this->current_scan_.scan_components_specification_parameters_.at((unsigned char)component_number).second;

        //We decode the DC component.
        decoded_dc = this->DecodeBaseline(file_content, index, this->dc_huffman_tables_.at(dc_table_index));
        diff = this->ReceiveBaseline(decoded_dc);
        diff = this->ExtendedBaseline(diff, decoded_dc);

        new_block.at<uchar>(0, 0, 0) = diff;

        // We decode the ac components.
        this->DecodeACCoefficients(file_content, index, &new_block, this->ac_huffman_tables_.at(ac_table_index));

        if (component_number == this->current_frame_header_.number_of_image_component)
        {
          component_number = 1;
        }
        else
        {
          component_number += 1;
        }
      }
      n += 1;
    }
  }
}

/**
 * \fn int JPEGDecoder::ReceiveBaseline(unsigned char number_of_bits)
 * \brief This function returns number_of_bits bits from the stream. The information is encoded as signed int.
 * 
 * \param[in] numebr_of_bits The number of bits to retrieve from the stream.
 */
int JPEGDecoder::ReceiveBaseline(unsigned char number_of_bits)
{
  int value = 0, i = 0;

  while (i != number_of_bits)
  {
    i += i;
    value = (value << 1) + this->NextBit();
  }

  return value;
}

/**
 * \fn int JPEGDecoder::ExtendedBaseline(unsigned char diff, unsigned char decoded_dc)
 * \brief The extend function returns the value that was encoded taking into account the range of encoding.
 * 
 * \param[in] difference The semi-encoded difference (depending on its sign).
 * \param[in] ssss The range of the encoded difference;
 */
int JPEGDecoder::ExtendedBaseline(int diff, unsigned char ssss)
{
  int value = 1;

  value = value << (ssss - 1);
  if (diff < value)
  {
    value = (-1 << ssss) + 1;
    diff += value;
  }
  return value;
}

void JPEGDecoder::DecodeRestartIntervalProgressive(unsigned char this->current_file_content_, int this->current_index_)
{
  this->ResetDecoderProgressive();

  while (!this->IsMarker(file_content, this->current_index_))
  {
    this->DecodeMCUProgressive(file_content, index);
  }
}

void JPEGDecoder::DecodeMCUProgressive(unsigned char this->current_file_content_, int this->current_index_)
{
  int n = 0;

  while (n < this->data_unit_per_mcu_)
  {
    // ADD the decoding function.
  }
}

void JPEGDecoder::InterpretFrameHeader(unsigned char this->current_file_content_, int this->current_index_, unsigned char encoding_process_type)
{
  // unsigned int header_length;
  unsigned char number_of_image_component_in_frame, mask_h = 240, mask_v = 15, c, horizontal_sampling_factor, vertical_sampling_factor, quantization_table_destination_selector;

  this->current_frame_header_.encoding_process_type_ = encoding_process_type;

  // To use to check header size.
  // header_length = int(file_content[this->current_index_] << 8 | file_content[this->current_index_ + 1]);
  this->current_frame_header_.sample_precision_ = file_content[this->current_index_ + 2];
  this->current_frame_header_.number_of_lines_ = int(file_content[this->current_index_ + 3] << 8 | file_content[this->current_index_ + 4]);
  this->current_frame_header_.number_of_samples_per_line_ = int(file_content[this->current_index_ + 5] << 8 | file_content[this->current_index_ + 6]);
  this->current_frame_header_.number_of_image_component = file_content[this->current_index_ + 7];
  this->current_index_ += 8;

  for (size_t i = 0; i < this->current_frame_header_.number_of_image_component; i++)
  {
    std::vector<unsigned char> components_vector;
    c = file_content[this->current_index_];
    horizontal_sampling_factor = (file_content[this->current_index_ + 1] & mask_h) >> 4;
    vertical_sampling_factor = file_content[this->current_index_ + 1] & mask_v;
    quantization_table_destination_selector = file_content[this->current_index_ + 2];
    components_vector.push_back(horizontal_sampling_factor);
    components_vector.push_back(vertical_sampling_factor);
    components_vector.push_back(quantization_table_destination_selector);
    this->current_frame_header_.component_signification_parameters_.insert(std::pair<unsigned char, std::vector<unsigned char>>(c, components_vector));
    this->current_index_ += 3;
  }
}

void JPEGDecoder::ParseHuffmanTableSpecification(unsigned char this->current_file_content_, int this->current_index_)
{
  unsigned int table_definition_length, counter, temp, sum_code_length = 0;
  unsigned char table_destination_identifier, table_class, mask_high = 240, mask_low = 15, number_of_huffman_codes;

  HuffmanTable table_being_parsed = {};

  table_definition_length = int(file_content[this->current_index_] << 8 | file_content[this->current_index_ + 1]);
  table_definition_length -= 2;
  this->current_index_ += 2;

  while (table_definition_length > 0)
  {
    sum_code_length = 0;
    table_class = (file_content[this->current_index_] & mask_high) >> 4;
    table_destination_identifier = file_content[this->current_index_] & mask_low;

    this->current_index_ += 1;

    for (size_t i = 0; i < 16; i++)
    {
      number_of_huffman_codes = file_content[this->current_index_];
      table_being_parsed.bits.push_back(number_of_huffman_codes);
      sum_code_length += number_of_huffman_codes;
      this->current_index_ += 1;
    }

    for (size_t i = 0; i < 16; i++)
    {
      counter = table_being_parsed.bits.at(i);
      table_being_parsed.huffvals.push_back(std::vector<unsigned char>());

      for (size_t j = 0; j < counter; j++)
      {
        table_being_parsed.huffvals.at(i).push_back(file_content[this->current_index_]);
        this->current_index_ += 1;
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
  table_being_parsed.huffcode = this->GenerateCodeTable(table_being_parsed.huffsize);
  this->DecoderTables(&table_being_parsed);

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

void JPEGDecoder::ParseScanHeader(unsigned char this->current_file_content_, int this->current_index_)
{
  // unsigned int scan_header_length;
  unsigned char image_component_in_scan, scan_component_selector, dc_table_destination_selector, ac_table_destination_selector, mask_high = 240, mask_low = 15;
  ScanHeader current_scan = {};

  // To do Check on the length of the header.
  // scan_header_length = int(file_content[this->current_index_] << 8 | file_content[this->current_index_ + 1]);
  this->current_index_ += 2;
  image_component_in_scan = file_content[this->current_index_];
  this->current_index_ += 1;

  this->current_scan_.number_of_image_components_ = image_component_in_scan;

  for (size_t i = 0; i < image_component_in_scan; i++)
  {
    scan_component_selector = file_content[this->current_index_];
    dc_table_destination_selector = (file_content[this->current_index_ + 1] & mask_high) >> 4;
    ac_table_destination_selector = file_content[this->current_index_ + 1] & mask_low;
    this->current_index_ += 2;
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

  this->current_scan_.start_of_spectral_selection_ = file_content[this->current_index_];
  this->current_scan_.end_of_spectral_selection_ = file_content[this->current_index_ + 1];
  this->current_scan_.approximation_high_bit_ = (file_content[this->current_index_ + 2] & mask_high) >> 4;
  this->current_scan_.approximation_low_bit_ = file_content[this->current_index_ + 2] & mask_low;
}

void JPEGDecoder::ParseComment(unsigned char this->current_file_content_, int this->current_index_)
{
  unsigned int comment_length;

  comment_length = (unsigned int)(file_content[this->current_index_] << 8 | file_content[this->current_index_ + 1]);

  std::cout << "Comment in the file : ";
  for (size_t i = 0; i < comment_length - 2; i++)
  {
    std::cout << file_content[this->current_index_ + 2 + i];
  }
  std::cout << std::endl;

  this->current_index_ += comment_length;
}

unsigned char JPEGDecoder::NextBit()
{
  unsigned char current_byte = file_content[this->current_index_], bit;

  if (this->next_bit_count_ == 0)
  {
    this->next_bit_count_ = 8;
    if (current_byte == 0xFF)
    {
      if (!(file_content[this->current_index_ + 1] == 0x00))
      {
        if (file_content[this->current_index_ + 1] == DEFINE_NUMBER_OF_LINE)
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
          this->current_index_ += 1;
        }
      }
    }
  }

  bit = current_byte << (8 - this->next_bit_count_);
  bit = bit >> 7;
  this->next_bit_count_ -= 1;

  if (this->next_bit_count_ == 0)
  {
    this->current_index_ += 1;
  }
  return bit;
}

/**
 * \fn unsigned char JPEGDecoder::DecodeBaseline(unsigned char this->current_file_content_, int this->current_index_, HuffmanTable used_table)
 * \brief This function reads the huffman code in the datastream and returns the huffval associated with the huffman code.
 * 
 * \param[in] file_content A pointer to the bitstream.
 * \param[in, out] index A pointer to the position of the cursor in the stream.
 * \param[in] used_table The table used for the decoding of the huffman code.
 */
unsigned char JPEGDecoder::DecodeBaseline(unsigned char this->current_file_content_, int this->current_index_, HuffmanTable used_table)
{
  int i = 1, j;
  char code;

  code = this->NextBit(file_content, index);
  while (code > used_table.max_code.at(i))
  {
    i += 1;
    code = (code << 1) + this->NextBit(file_content, index);
  }

  j = used_table.val_pointer.at(i - 1);
  j = j + code - used_table.max_code.at(i - 1);
  return used_table.huffvals.at(i - 1).at(j);
}

/**
 * \fn std::vector<unsigned char> JPEGDecoder::GenerateSizeTable(std::vector<unsigned char> bits)
 * \brief Function to get the huffsize table for a huffman tree. The huffsize table contains all the size of the codes in the huffman table, i.e if we have three codes of length 3 and to of length 4, the huffsize will be [3,3,3,4,4].
 * 
 * \param[in] bits This is the array containing the number of code of each length, starting from 1.
 */
std::vector<unsigned char> JPEGDecoder::GenerateSizeTable(std::vector<unsigned char> bits)
{
  unsigned char k = 0, i = 1, j = 1;
  std::vector<unsigned char> huffsize;

  huffsize.push_back(0);

  while (i < 16)
  {
    if (j > bits.at(i - 1))
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
  huffsize.push_back(0);
  this->last_k_ = k;
  return huffsize;
}

/**
 * \fn void JPEGDecoder::GenerateCodeTable(HuffmanTable *table_to_fill)
 * \brief Function to generate the code table of the huffman tree. The code are generated using the huffsize table generated from bits.
 * 
 * \param[in] huffsize The array containing the list of the sizes of all the codes.
 */
std::vector<unsigned short> JPEGDecoder::GenerateCodeTable(std::vector<unsigned char> huffsize)
{
  int k = 0, code = 0;
  unsigned char si = huffsize.at(0);
  std::vector<unsigned short> huffcode(huffsize.size(), 0);

  while (true)
  {
    do
    {
      huffcode.at(k) = code;
      code += 1;
      k += 1;
    } while (huffsize.at(k) == si);

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

/**
 * \fn void JPEGDecoder::DecoderTables(HuffmanTable *table_processed)
 * \brief For a given huffman table, this function creates the table MAXCODE, MINCODE and VALPTR, which respectively hold the maximum code in value for a given size, the minimum code for a given size and the pointer to these values in the other tables (huffval, huffsize). The three table generated hold 17 values to have the index starting from 1 (match the norm representation).
 * 
 * \param[in, out] table_processed Pointer to the huffman table being currently processed.
 */
void JPEGDecoder::DecoderTables(HuffmanTable *table_processed)
{
  int i = 0, j = 0;
  std::vector<short> max_code(17, 0), min_code(17, 0), val_ptr(17, 0);

  while (true)
  {
    i += 1;
    if (i > 16)
    {
      return;
    }
    if (table_processed->bits[i - 1])
    {
      max_code.at(i) = -1;
    }
    else
    {
      val_ptr.at(i) = j;
      min_code.at(i) = table_processed->huffcode.at(j);
      j = j + table_processed->bits[i - 1] - 1;
      max_code.at(i) = table_processed->huffcode.at(j);
      j += 1;
    }
  }
}

void JPEGDecoder::ResetDecoderBaseline()
{
}

bool JPEGDecoder::IsMarker(unsigned char this->current_file_content_, int index)
{

  if (file_content[index] == 0xFF)
  {
    if (file_content[index + 1] == 0xFF)
    {
      return false;
    }
    return true;
  }
  else
  {
    return false;
  }
}

/**
 * \fn void JPEGDecoder::DecodeACCoefficients(unsigned char this->current_file_content_, int this->current_index_, cv::Mat *new_block)
 * \brief Decode the AC coefficients for the baseline procedure and store them in the new_block.
 * 
 * \param[in] file_content Pointer to the file being parsed.
 * \param[in, out] index, Pointer to the position of the cursor in the file being parsed.
 * \param[in, out] new_block Pointer to the block of data to contain the decoded values.
 */
void JPEGDecoder::DecodeACCoefficients(unsigned char this->current_file_content_, int this->current_index_, cv::Mat *new_block, HuffmanTable used_table)
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
      ZZ.at(k) = this->DecodeZZ(ssss);
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
unsigned char JPEGDecoder::DecodeZZ(unsigned char ssss)
{
  unsigned char return_value;
  return_value = this->ReceiveBaseline(ssss);
  return_value = this->ExtendedBaseline(return_value, ssss);
  return return_value;
}

/**
 * \fn unsigned char *JPEGDecoder::GetMarker()
 * \brief Returns the current marker if any, else throws a runtime error.
 */
unsigned char *JPEGDecoder::GetMarker()
{
  std::stringstream error;
  unsigned char *marker;

  if (this->current_file_content_[this->current_index_] == 0xFF)
  {
    marker = new unsigned char[1];
    std::memcpy(marker, &(this->current_file_content_[this->current_index_ + 1]), 1);
    this->current_index_ = this->current_index_ + 2;
    return marker;
  }
  else
  {
    std::cout << this->current_index_ << std::endl;
    error << "Error while reading marker, 0xFF expected, but " << std::hex
          << (int)this->current_file_content_[this->current_index_] << " found.";
    throw std::runtime_error(error.str());
  }
}