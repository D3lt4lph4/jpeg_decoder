# JPEG Baseline Decoder

This project aims to create a "level" depending decoder that can decode an image and stop at a given point in the decoding.
The only supported process for decoding is the baseline.

## How does the compression works.

The workflow of the compression/decompression is the following one (main steps), the image was taken from wikipedia :

![compression/decompression workflow](https://raw.githubusercontent.com/D3lt4lph4/jpeg_encoder_decoder/master/images/compression_JPEG.png?token=AXSrihw6StMXldgUNoZ5d55DTkqKOXrGks5bdYj-wA%3D%3D "JPEG workflow")


The following parts will describe the process for each steps in the image above.

## Color space

The images are mostly handled in RGB (or BGR) by the image processing framework/library since it is easier to handle for a human. But for storing the data, the YCbCr space is better, the human eyes does not have the same sensibility to the different component thus, some under-sampling can be carried on the data in this space.

## UnderSampling

## DCT

## Quantification 

## Encoding