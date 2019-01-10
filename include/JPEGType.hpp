#ifndef __JPEG_TYPE__
#define __JPEG_TYPE__

#include <map>
#include <vector>

const unsigned char ZZ_order[] = {
    0,  1,  5,  6,  14, 15, 27, 28, 2,  4,  7,  13, 16, 26, 29, 42,
    3,  8,  12, 17, 25, 30, 41, 43, 9,  11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54, 20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61, 35, 36, 48, 49, 57, 58, 62, 63};
const unsigned char FRAME_TYPE_BASELINE_DTC = 0xB0;
const unsigned char FRAME_TYPE_PROGRESSIVE = 0xB2;

const unsigned char START_OF_FRAME_BASELINE = 0xC0;
const unsigned char START_OF_FRAME_PROGRESSIVE = 0xC2;
const unsigned char DEFINE_HUFFMAN_TABLE = 0xC4;

const unsigned char START_OF_IMAGE = 0xD8;
const unsigned char END_OF_IMAGE = 0xD9;
const unsigned char START_OF_SCAN = 0xDA;
const unsigned char DEFINE_QUANTIZATION_TABLE = 0xDB;
const unsigned char DEFINE_NUMBER_OF_LINE = 0xDC;
const unsigned char DEFINE_RESTART_INTERVAL = 0xDD;
const unsigned char DEFINE_HIERARCHICAL_PROGRESSION = 0xDE;
const unsigned char EXPAND_REFERENCE_COMPONENTS = 0xDF;

const unsigned char COMMENT = 0xFE;

// Define all the APP blocks
const unsigned char APPO = 0xE0;
const unsigned char APP1 = 0xE1;
const unsigned char APP2 = 0xE2;
const unsigned char APP3 = 0xE3;
const unsigned char APP4 = 0xE4;
const unsigned char APP5 = 0xE5;
const unsigned char APP6 = 0xE6;
const unsigned char APP7 = 0xE7;
const unsigned char APP8 = 0xE8;
const unsigned char APP9 = 0xE9;
const unsigned char APP10 = 0xEA;
const unsigned char APP11 = 0xEB;
const unsigned char APP12 = 0xEC;
const unsigned char APP13 = 0xED;
const unsigned char APP14 = 0xEE;
const unsigned char APP15 = 0xEF;

const unsigned char JFIF[] = {0x4a, 0x46, 0x49, 0x46, 0x00};

enum ImageType {blackandwhite, rgb, rgba};

/**
 * \struct QuantizationTable
 * \brief The representation of the quantization table defined in the JPEG norm.
 *
 * \var QuantizationTable::pq_
 * Quantization table element precision. Specifies
 * the precision of the Qk values. Value 0 indicates 8-bit Qk values; value 1
 * indicates 16-bit Qk values. Pq shall be zero for 8 bit sample precision P.
 *
 * \var QuantizationTable::qks_
 * Vector containing all the values for the
 * quantization table, should be of size 64.
 */
struct QuantizationTable {
  unsigned char pq_;
  std::vector<unsigned short> qks_;
};

/**
 * \struct JFIFHeader
 * \brief The representation of a JFIF Header. Contains all the usefull
 * variables found after the JFIF marker.
 *
 * \var JFIFHeader::current_version_
 * The version of the revision of the JFIF norm used when encoding the file.
 *
 * \var JFIFHeader::current_unit_
 * The unit to use for the pixel density, a value of 0 means that pixel_density
 * values are  aspect ratio, a value of 1 means dots per inch and a value of 2
 * means dots per centimeters.
 *
 * \var JFIFHeader::horizontal_pixel_density_
 * The horizontal pixel density for the image, the unit used is defined by
 * current_unit_.
 *
 * \var JFIFHeader::vertical_pixel_density_
 * The vertical pixel density for the image, the unit used is defined by
 * current_unit_.
 *
 * \var JFIFHeader::thumbnail_horizontal_pixel_count_
 * The horizontal pixel count for the thumbnail if define.
 *
 * \var JFIFHeader::thumbnail_vertical_pixel_count_
 * The vertical pixel count for the thumbnail if define.
 */
struct JFIFHeader {
  unsigned short current_version_;
  unsigned char current_unit_, horizontal_pixel_density_,
      vertical_pixel_density_, thumbnail_horizontal_pixel_count_,
      thumbnail_vertical_pixel_count_;
};

/**
 * \struct FrameHeader
 * \brief This structure regroup the informations about the frame header. The
 * frame header is here to specify the informations about the source image.
 *
 * \var FrameHeader::encoding_process_type_
 * The decoding process to use. Only the baseline is supported for now.
 *
 * \var FrameHeader::sample_precision_
 * Specifies the precision in bits for the samples of the components in the
 * frame.
 *
 * \var FrameHeader::number_of_lines_
 * Specifies the maximum number of lines in the source image. See the norm for
 * more details.
 *
 * \var FrameHeader::number_of_samples_per_line_
 * Specifies the maximum number of samples line in the source image. See the
 * norm for more details.
 *
 * \var FrameHeader::number_of_component_
 * The number of components in the source image.
 *
 * \var FrameHeader::component_parameters_
 * A map containing the parameters for each components. The key is the
 * identifier used in the coding process the values are H the horizontal
 * sampling factor, V the vertical sampling factor and Tq the quantization table
 * destination selector (i.e which quantization table use).
 */
struct FrameHeader {
  unsigned char encoding_process_type_, sample_precision_;
  unsigned int number_of_lines_, number_of_samples_per_line_,
      number_of_component_;
  std::map<unsigned char, std::vector<unsigned char>> component_parameters_;
};

/**
 * \struct ScanHeader
 * \brief Contains the information about a scan.
 *
 * \var ScanHeader::number_of_component_s_
 * The number of image components in the scan.
 *
 * \var ScanHeader::start_of_spectral_selection_
 * Specifies the first DCT coefficient in each block in zig-zag order which
 * shall be coded in the scan.
 *
 * \var ScanHeader::end_of_spectral_selection_
 * Specifies the last DCT coefficient in each block in zig-zag order which shall
 * be coded in the scan.
 *
 * \var ScanHeader::approximation_high_bit_
 * Specifies the point transform used in the preceding scan.
 *
 * \var ScanHeader::approximation_low_bit_
 * Specifies the point transform (should be 0 in baseline mode).
 *
 * \var ScanHeader::components_parameters_
 * A map containing the parameters for each components. The key is the
 * identifier used in the coding process, the values are Td the DC entropy table
 * to be used for the component and Ta the AC entropy table to be used for the
 * component.
 */
struct ScanHeader {
  unsigned char number_of_component_s_, start_of_spectral_selection_,
      end_of_spectral_selection_, approximation_high_bit_,
      approximation_low_bit_;
  std::map<unsigned char, std::pair<unsigned char, unsigned char>>
      components_parameters_;
};

/**
 * \struct HuffmanTable
 * \brief Contains the representation of a Huffman table.
 *
 * \var HuffmanTable::table_class_
 * Either 0 or 1 respectively for DC and AC tables.
 *
 * \var HuffmanTable::last_k_
 * Value used in the generate size table algorithm.
 *
 * \var HuffmanTable::bits
 * A list containing the number of code of each size (from 1 to 16), i.e [0,2,4,0,0,...] would be two code of size 2 and 4 of size 4.
 *
 * \var HuffmanTable::huffsize
 * A list containing each code length i.e if we have three codes of length 3 and to of length 4, the huffsize will be [3,3,3,4,4,0]. The 0 terminates the list.
 *
 * \var HuffmanTable::val_pointer
 * A list containing the index of the min_code value in the huffcode list.
 *
 * \var HuffmanTable::huffvals
 * The symbol value associated with each code length. This is the number of bits to read in the stream.
 *
 * \var HuffmanTable::huffcode
 * A list containing the code for each of the length, i.e with huffsize of [3,3,3,4,4,0] we would have something like (not actual values) [100, 101, 111, 1001, 1010].
 *
 * \var HuffmanTable::min_code
 * The smallest code for a given size (given 111, 101, 100, min_code[2] == 100 )
 *
 * \var HuffmanTable::max_code
 * The biggest code for a given size (given 111, 101, 100, max_code[2] == 111 )
 */
struct HuffmanTable {
  unsigned char table_class_, last_k_;
  std::vector<unsigned char> bits, huffsize, val_pointer, huffvals;
  std::vector<unsigned short> huffcode;
  std::vector<int> min_code, max_code;
};
#endif