""" This script aims to parse a whole dataset to see if the processing of any image is bugged.
"""

import argparse
from os.path import isfile
import glob

from tqdm import tqdm

from jpegdecoder.decoder import JPEGDecoder

parser = argparse.ArgumentParser()

parser.add_argument("directory", help="The head directory to be processed. If the recursive option is specified, only the leaves of the tree will be checked for images.")
parser.add_argument("-r", "--recursive", help="Tells if the script should process the subfolder in the head folder. Only the leaves will be checked for images if the option is specified.")
parser.add_argument("-o", "--output", help="The name on the file that will contain all the names of the bugged images.", default="bugged.txt")

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

for image in tqdm(res):
    try:
        if isfile(image):
            _ = decoder.decode_file(image, 4)
    except Exception as e:
        with open(args.output, "a+") as file_bugged_images:
            file_bugged_images.write("{}\n".format(image))
