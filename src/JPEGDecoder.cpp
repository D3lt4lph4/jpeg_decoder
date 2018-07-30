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
        case DEFINE_QUANTIZATION_TABLE:
          this->ParseQuantizationTable(file_content, &current_index);
          std::cout << "quantization table parsed" << std::endl;
          break;
        case START_OF_FRAME_PROGRESSIVE:
          this->ParseFrameHeader(file_content, &current_index);
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
  QuantizationTable *table_being_parsed, insert_res;

  lq = int(file_content[*index] << 8 | file_content[*index + 1]);
  lq -= 2;
  *index += 2;
  while (lq > 0)
  {
    table_being_parsed = new QuantizationTable();
    pq = file_content[*index];
    tq = (pq & mask_tq) > 4;
    pq = pq & mask_pq;

    table_being_parsed->pq_ = pq;

    *index += 1;
    if (pq == 1)
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[*index] << 8 | file_content[*index + 1]);
        *index += 2;
        table_being_parsed->qks_.push_back(qk);
      }
    }
    else
    {
      for (size_t i = 0; i < 64; i++)
      {
        qk = int(file_content[*index]);
        *index += 1;
        table_being_parsed->qks_.push_back(qk);
      }
    }
    if (!this->quantization_tables_
             .insert(std::pair<unsigned char, QuantizationTable *>(
                 tq, table_being_parsed))
             .second)
    {
      this->quantization_tables_.erase(tq);
      this->quantization_tables_.insert(
          std::pair<unsigned char, QuantizationTable *>(tq,
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

void JPEGDecoder::ParseFrameHeader(unsigned char *file_content, int *index)
{
  unsigned int lf, y, x;
  unsigned char p, nf, mask_h = 240, mask_v = 15, c, h, v, tq;

  lf = int(file_content[*index] << 8 | file_content[*index + 1]);
  p = file_content[*index + 2];
  lf = int(file_content[*index + 3] << 8 | file_content[*index + 4]);
  lf = int(file_content[*index + 5] << 8 | file_content[*index + 6]);
  nf = file_content[*index + 7];
  *index += 8;

  for (size_t i = 0; i < nf; i++)
  {
    c = file_content[*index];
    h = (file_content[*index + 1] & mask_h) >> 4;
    v = file_content[*index + 1] & mask_v;
    tq = file_content[*index + 2];
    *index += 3;
  }
}

void JPEGDecoder::ParseHuffmanTableSpecification(unsigned char *file_content, int *index)
{
  unsigned int lh, counter, temp, sum_l = 0;
  unsigned char tc, th, mask_tc = 240, mask_th = 15, l;
  std::vector<unsigned char> ls, vs;

  lh = int(file_content[*index] << 8 | file_content[*index + 1]);
  lh -= 2;
  *index += 2;

  while (lh > 0)
  {
    sum_l = 0;
    tc = (file_content[*index + 1] & mask_tc) >> 4;
    th = file_content[*index + 1] & mask_th;
    *index += 1;
    for (size_t i = 0; i < 16; i++)
    {
      l = file_content[*index];
      ls.push_back(l);
      sum_l += l;
      *index += 1;
    }

    for (size_t i = 0; i < 16; i++)
    {
      counter = ls.at(i);

      for (size_t j = 0; j < counter; j++)
      {
        vs.push_back(file_content[*index]);
        *index += 1;
      }
    }
    lh = lh - 17 - sum_l;
  }
}

void JPEGDecoder::ParseScanHeader(unsigned char *file_content, int *index)
{
  unsigned int ls;
  unsigned char ns, cs, td, ta, ss, se, ah, al, mask_td = 240, mask_ta = 15;

  ls = int(file_content[*index] << 8 | file_content[*index + 1]);
  *index += 2;
  ns = file_content[*index];
  *index += 1;

  for (size_t i = 0; i < ns; i++)
  {
    cs = file_content[*index];
    td = (file_content[*index + 1] & mask_td) >> 4;
    ta = file_content[*index + 1] & mask_ta;
    *index += 2;
  }

  ss = file_content[*index];
  se = file_content[*index + 1];
  ah = (file_content[*index + 2] & mask_td) >> 4;
  al = file_content[*index + 2] & mask_ta;
}