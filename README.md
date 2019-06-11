# JPEG Baseline Decoder

This project aims to create a "level" depending decoder that can decode an image and stop at a given point in the decoding pipeline.
The only supported process for decoding is the baseline. Only the versions starting from 2.2.0 should be used, anything below should work but is not documented and incomplete. Starting from 2.2.0 the pipeline to install and use the decoder will remain the same. The project was done in c++ and can be use in python thanks to [pybind11](https://github.com/pybind/pybind11)

This project was first done as part of my thesis, the aim was to access the DCT coefficient from the JPEG compression pipeline. As of today, this is more like a side project and I recommend using [this decoder](https://github.com/uber-research/jpeg2dct) if someone was to access the coefficient, it is faster, lighter and probably simpler to use.

## How does the compression works ?

The workflow of the compression/decompression is the following one (main steps) :

![compression/decompression workflow](https://github.com/D3lt4lph4/jpeg_decoder/blob/master/images/JPEG_pipeline_full.png "JPEG workflow")

The first two steps are the YCbCr transform and the sub-sampling, the space transformation is done to be able to sub-sample without loosing too much information (the human eye is less sensitive to variation in the CbCr components). Then the DCT and quantization transformations are performed to further compress the data. The high frequencies are "set to 0" through the quantization to take advantage of the compression RLE/Huffman.

## Simple install in a python

For all the steps below, I will assume you have set up and activated a virtualenv. First let's show the simplest installation:

```bash
python setup.py install
pip install numpy
# or python setup.py sdist
```

To test in python

```python
import jpegdecoder
import numpy as np

decoder = jpegdecoder.decoder.JPEGDecoder()
# The 2 is the decoding level
image = decoder.decode_file("data/chat_420.jpg", 2)

rows, cols = image.get_component_shape(0)[0:2]
y = np.reshape(image.get_data(0), image.get_component_shape(0))
Cb = np.reshape(image.get_data(1), image.get_component_shape(1))
Cr = np.reshape(image.get_data(2), image.get_component_shape(2))
```

The decode_file function takes the path to the file and the decoding level as input. The decoding level is as follow:

- 1: the dequantization step is not performed
- 2: the IDCT step is not performed
- 3: the YCbCr to RGB step is not performed
- 4: returns the image

## Compile the project

This is for compilation of the c++ code. Starting from version 2.2.0 the project has two different mode for compilation DEBUG and RELEASE. Depending on the mode chosen, the feature are not exactly the same.

Debug:

- All the test are compiled and can be run to check if everything is working
- The main file contains opencv support to display the decompressed images
- No python wrapper is output for now

Release:

- The test are not compiled
- No opencv support (the images cannot be displayed)
- The python wrappers are created as well as the package to distribute them

To compiling of the project is almost the same in both case :

```bash
# Compiling in debug mode
mkdir debug
cd debug
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
make

# Will fully decode the image and show it
./jpeg_decoder -f ../data/chat_420.jpg -s

# Running the tests
./test/test_huffman
./test/test_parsing
./test/test_utility

# Compiling in release mode
git submodule init
git submodule update
mkdir release
cd release
cmake ..
make
```

## Documentation

The documentation is not yet finished, to have more details for the python available function go check the _python.cpp files in the src folder.

## Known Bugs

### OpenCV

OpenCV seems to look for cuda while compiling, if you have multiple version of cuda installed on the computer, cmake might stop at the first one encountered and not look further even if this not the right version. To remedy to this problem, rename the wrong cuda-X-X directories in /usr/local to something else (ugly fix but working).

## Version description

### v1.0.0

Features:

- Can decompress at various steps of the JPEG deocing process
- The data can be saved into a file to later use with other language
- python script provided to read the saved data
- Logging information available

Limitations:

- The implementation of the IDCT is naive and thus slow
- The data is not saved under the subsampled form
- The 8x8 DCT blocks are duplicated instead of being correctly upsampled
- Usage of OpenCV data structure for the "images"

### v1.1.0

Features:

- Implementation of the fast IDCT to improve the decoding speed.
- Setting the output name of the files depending on the step of decompression
  - .qhjpg for de huffman and dequantized data
  - .iqhjpg for de huffman, dequantized data and IDCT
  - .riqhjpg for .jpg (? should be changed)

Limitations:

- Usage of OpenCV data structure for the "images"
- The data is not saved under the subsampled form
- The 8x8 DCT blocks are duplicated instead of being correctly upsampled

### v1.2.0

Features:

- Removing dependency to OpenCV for the decoder
- OpenCV included in the main file only in debug mode to display the results
- Adding documentation

Limitations:

- The data is not saved under the subsampled form
- The 8x8 DCT blocks are duplicated instead of being correctly upsampled

### v2.0.0

Features:

- The image is correctly returned in the subsampled form
- Adding a new type, JPEGImage to store the data

Limitations:

Warnings:

- The data is no longer saved in the same way, incompatibility between the new and old python decoding functions

### v2.1.0

Features:

- Improving the speed of the decoder through the usage of reference and (smart)pointers

### v2.2.0

Features:

- Adding wrapper for python usage of the library
    - Numpy array are used as output for the data of the images
    - Documentation for the Python functions available
    - Testing of the python side of the library
- Documentation for each of the tagged version was added on the ReadTheDocs server.
- Boost is removed from compilation if not in debug mode
- Improved tests for the library on the c++ side
- Adding files to calculate the time the decoder takes to decode images at different levels
