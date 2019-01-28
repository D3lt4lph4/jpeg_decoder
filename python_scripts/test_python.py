from jpegdecoder import decoder
import sys
dec = decoder.JPEGDecoder()
img = dec.decode_file("/data/datasets/imagenet/ILSVRC_2012_224_224/validation/n01440764/ILSVRC2012_val_00000293.JPEG", 2)
img.get_data(0)
del dec
img.get_data(0)
