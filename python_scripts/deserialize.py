import cv2 as cv
import numpy as np

f = open("output.txt", "rb")

try:
    rows = int.from_bytes(f.read(4), "little")
    cols = int.from_bytes(f.read(4), "little")
    mat_type = int.from_bytes(f.read(4), "little")
    channel = int.from_bytes(f.read(4), "little")
    
    res = np.zeros((rows, cols, channel), dtype="int32")
    for row in range(rows): 
        for col in range(cols):
            for chan in range(channel):
                res[row, col, chan] = int.from_bytes(f.read(4), "little")
finally:
    f.close()

cv.imshow('', res)
cv.waitKey()
