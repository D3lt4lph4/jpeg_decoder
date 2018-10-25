import cv2 as cv2
import numpy as np

f = open("datachat_420.qhjpg", "rb")

try:
    # First we read the number of components
    components_number = int.from_bytes(f.read(4), "little")
    print(components_number)

    shapes = np.zeros((components_number, 2), dtype="int32")
    
    # Then we read the shapes into a numpy array for each of the components. We store them because the shapes might be different for each components if we are not in RGB space.
    for i in range(components_number):
        shapes[i, 0] = int.from_bytes(f.read(4), "little")
        shapes[i, 1] = int.from_bytes(f.read(4), "little")

    real_shape = np.zeros(components_number, dtype="int32")
    # Then we read the real image size
    for i in range(components_number):
        real_shape[i] = int.from_bytes(f.read(4), "little")

    # Finally we read the data from the stream

    res = []
    for i in range(components_number):
        print("{}, {}".format(shapes[i, 0], shapes[i, 1]))
        data = np.zeros((shapes[i, 0], shapes[i, 1]), dtype="int32")
        for row in range(shapes[i, 0]):
            for col in range(shapes[i, 1]):
                data[row, col] = int.from_bytes(f.read(4), "little")
        res.append(data)

    print(res[0].shape)
    print(res[1].shape)
    print(res[2].shape)

    # res = np.array(res, dtype = np.uint8)
    # res = np.swapaxes(res, 0, 2)
    # res = np.swapaxes(res, 0, 1)
    # cv2.imshow('Color image', res)
    # cv2.waitKey(0)

finally:
    f.close()
