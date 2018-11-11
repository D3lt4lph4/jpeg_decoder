# JPEG Baseline Decoder

This project aims to create a "level" depending decoder that can decode an image and stop at a given point in the decoding.
The only supported process for decoding is the baseline.

## How does the compression works.

The workflow of the compression/decompression is the following one (main steps), the image was taken from wikipedia :

![compression/decompression workflow](https://raw.githubusercontent.com/D3lt4lph4/jpeg_encoder_decoder/master/images/compression_JPEG.png?token=AXSrihw6StMXldgUNoZ5d55DTkqKOXrGks5bdYj-wA%3D%3D "JPEG workflow")

# Version tag description

## v1.0.0

First working version of the decoder:

- The implementation of the IDCT is naive and thus slow
- The 8x8 DCT blocks are duplicated instead of being correctly upsampled

## v1.1.0

Features added:

- Implementation of the fast IDCT to improve the decoding speed.

## v1.2.0

Removing OpenCV from the project if not in DEBUG mode.

## v2.0.0



## v2.1.0

Removing OpenCV from the project if not in DEBUG mode.
