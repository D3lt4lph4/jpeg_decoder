""" This script aims to parses x images and output a file with all the matrix with the DCT coefficient and matching YCbCr IDCT results.
"""

import argparse
from os.path import isfile
import glob
import numpy as np

from tqdm import tqdm

from jpegdecoder.decoder import JPEGDecoder

parser = argparse.ArgumentParser()

parser.add_argument("directory", help="The head directory to be processed. If the recursive option is specified, only the leaves of the tree will be checked for images.")
parser.add_argument("-r", "--recursive", help="Tells if the script should process the subfolder in the head folder. Only the leaves will be checked for images if the option is specified.")
parser.add_argument("-n", "--number", help="The number of images to process.", type=int , default=50)

args = parser.parse_args()

# Setting the recursive option
if args.recursive is not None:
    args.recursive == True
else:
    args.recursive == False

# find all the files in the dataset
dataset_path_template = "{}/**/*".format(args.directory)
res = glob.glob(dataset_path_template, recursive=args.recursive)

# Trying to decode all the images
decoder = JPEGDecoder()

# Creating the matrix to hold the results
Y_dct = []
Cb_dct = []
Cr_dct = []
Y = []
Cb = []
Cr = []
i = 0
for image in tqdm(res):
    if i == args.number:
        break
    if isfile(image):
        
        img_dct = decoder.decode_file(image, 2)
        if img_dct.get_real_shape()[2] == 1:
            continue
        rows, cols = img_dct.get_component_shape(0)
        i += 1
        Y_dct_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        Cb_dct_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        Cr_dct_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        img_dct_data_0 = img_dct.get_data(0)

        img_dct_data_1 = img_dct.get_data(1)
        img_dct_data_2 = img_dct.get_data(2)
        for row in range(rows):
            for col in range(cols):
                Y_dct_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_dct_data_0[row * cols + col]
                Cb_dct_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_dct_data_1[row * cols + col]
                Cr_dct_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_dct_data_2[row * cols + col]
        
        Y_dct.append(Y_dct_image)
        Cb_dct.append(Cb_dct_image)
        Cr_dct.append(Cr_dct_image)

        img = decoder.decode_file(image, 3)
        
        rows, cols = img.get_component_shape(0)

        Y_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        Cb_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        Cr_image = np.zeros((rows // 8 * cols // 8, 8, 8))
        img_data_0 = img.get_data(0)
        img_data_1 = img.get_data(1)
        img_data_2 = img.get_data(2)
        for row in range(rows):
            for col in range(cols):
                Y_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_data_0[row * cols + col]
                Cb_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_data_1[row * cols + col]
                Cr_image[col // 8 + row // 8 * (cols // 8), row % 8, col % 8] = img_data_2[row * cols + col]
        Y.append(Y_image)
        Cb.append(Cb_image)
        Cr.append(Cr_image)

Y_dct = np.array(Y_dct)
Cb_dct = np.array(Cb_dct)
Cr_dct = np.array(Cr_dct)
Y = np.array(Y)
Cb = np.array(Cb)
Cr = np.array(Cr)

np.savez("data.npz", Y_dct, Cb_dct, Cr_dct, Y, Cb, Cr)