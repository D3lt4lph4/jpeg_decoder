import numpy as np

f = open("output.txt", "rb")

try:
    rows = int.from_bytes(f.read(4), "little")
    cols = int.from_bytes(f.read(4), "little")
    channel = int.from_bytes(f.read(4), "little")

    res = np.zeros((rows, cols, channel), dtype="int32")

    for row in range(rows/8):
        for col in range(cols/8):
            for chan in range(channel):
                for row_cell in range(8):
                    for col_cell in range(8):
                        res[row*8+row_cell, col*8+col_cell, chan] = int.from_bytes(f.read(4), "little")
finally:
    f.close()
