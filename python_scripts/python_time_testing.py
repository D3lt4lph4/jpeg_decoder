import argparse
import time
from os import listdir
from os.path import isfile, join
import numpy as np

from tqdm import tqdm

import cv2

import jpegdecoder

parser = argparse.ArgumentParser()

parser.add_argument("directory", help="The directory containing all the images to use for the test.")
parser.add_argument("-o", "--output", help="The output file for the results of the test, if not specified the results are only printed in the console.")

args = parser.parse_args()

image_directory = args.directory

images =  [join(image_directory, f) for f in listdir(image_directory) if isfile(join(image_directory, f))]

results = {}

print("-------------------------Processing with OpenCV-------------------------")
print("Processing all the images with OpenCV (not stored in memory).")
start_time_opencv_not_kept = time.time()
for image in tqdm(images):
    _ = cv2.imread(image,0)
results["elapsed_time_opencv_not_kept"] = time.time() - start_time_opencv_not_kept

print("Processing all the images with OpenCV (stored in memory).")
images_results = []
start_time_opencv_kept = time.time()
for image in tqdm(images):
    images_results.append(cv2.imread(image,0))
    print(images_results[-1].shape)
results["elapsed_time_opencv_kept"] = time.time() - start_time_opencv_kept


print("-------------------------Processing with the jpegdecoder-------------------------")
print("Processing all the images with the JPEGDecoder, RGB (not stored in memory).")

decoder = jpegdecoder.decoder.JPEGDecoder()
start_time_jpegdecoder_not_kept = time.time()
for image in tqdm(images):
    jpg_image = decoder.decode_file(image, 4)
    img_real_size = jpg_image.get_real_shape()
    img = np.empty((img_real_size[0], img_real_size[1], img_real_size[2]))
    img[:,:,0] = np.resize(jpg_image.get_data(0), (img_real_size[0], img_real_size[1]))
    img[:,:,1] = np.resize(jpg_image.get_data(1), (img_real_size[0], img_real_size[1]))
    img[:,:,2] = np.resize(jpg_image.get_data(2), (img_real_size[0], img_real_size[1]))

results["elapsed_time_jpeg_rgb_not_kept"] = time.time() - start_time_jpegdecoder_not_kept

print("Processing all the images with the JPEGDecoder, RGB (stored in memory).")
decoder = jpegdecoder.decoder.JPEGDecoder()
images_results = []
start_time_jpegdecoder_kept = time.time()
for image in tqdm(images):
    jpg_image = decoder.decode_file(image, 4)
    img_real_size = jpg_image.get_real_shape()
    img = np.empty((img_real_size[0], img_real_size[1], img_real_size[2]))
    img[:,:,0] = np.resize(jpg_image.get_data(0), (img_real_size[0], img_real_size[1]))
    img[:,:,1] = np.resize(jpg_image.get_data(1), (img_real_size[0], img_real_size[1]))
    img[:,:,2] = np.resize(jpg_image.get_data(2), (img_real_size[0], img_real_size[1]))

    images_results.append(img)
results["elapsed_time_jpeg_rgb_kept"] = time.time() - start_time_jpegdecoder_kept

print("Processing all the images with the JPEGDecoder, DCT (not stored in memory).")
decoder = jpegdecoder.decoder.JPEGDecoder()
start_time_jpegdecoder_not_kept_dct = time.time()
for image in tqdm(images):
    jpg_image = decoder.decode_file(image, 2)
    img_real_size = jpg_image.get_real_shape()
    img = np.empty((img_real_size[0], img_real_size[1], img_real_size[2]))
    img[:,:,0] = np.resize(jpg_image.get_data(0), (img_real_size[0], img_real_size[1]))
    img[:,:,1] = np.resize(jpg_image.get_data(1), (img_real_size[0], img_real_size[1]))
    img[:,:,2] = np.resize(jpg_image.get_data(2), (img_real_size[0], img_real_size[1]))
results["elapsed_time_jpeg_dct_not_kept"] = time.time() - start_time_jpegdecoder_not_kept_dct

print("Processing all the images with the JPEGDecoder, DCT (stored in memory).")
decoder = jpegdecoder.decoder.JPEGDecoder()
images_results = []
start_time_jpegdecoder_kept_dct = time.time()
for image in tqdm(images):
    jpg_image = decoder.decode_file(image, 2)
    img_real_size = jpg_image.get_real_shape()
    img = np.empty((img_real_size[0], img_real_size[1], img_real_size[2]))
    img[:,:,0] = np.resize(jpg_image.get_data(0), (img_real_size[0], img_real_size[1]))
    img[:,:,1] = np.resize(jpg_image.get_data(1), (img_real_size[0], img_real_size[1]))
    img[:,:,2] = np.resize(jpg_image.get_data(2), (img_real_size[0], img_real_size[1]))
    images_results.append(img)
results["elapsed_time_jpeg_dct_kept"] = time.time() - start_time_jpegdecoder_kept_dct

print(results)
# print("The results are as follow: ")
# print("\t - OpenCV, RGB, images not stored in memory: {}".format(results["elapsed_time_opencv_not_kept"]))
# print("\t - OpenCV, RGB, images stored in memory: {}".format(results["elapsed_time_opencv_kept"]))
# print("\t - JPEGDecoder, RGB, images not stored in memory: {}".format(results["elapsed_time_jpeg_rgb_not_kept"]))
# print("\t - JPEGDecoder, RGB, images stored in memory: {}".format(results["elapsed_time_jpeg_rgb_kept"]))
# print("\t - JPEGDecoder, DCT, images not stored in memory: {}".format(results["elapsed_time_jpeg_dct_not_kept"]))
# print("\t - JPEGDecoder, DCT, images not stored in memory: {}".format(results["elapsed_time_jpeg_dct_kept"]))