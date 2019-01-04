import numpy as np

npzfile = np.load("data.npz")

print(npzfile.files)

Y_dct = npzfile['arr_0']
Cb_dct = npzfile['arr_1']
Cr_dct = npzfile['arr_2']
Y = npzfile['arr_3']
Cb = npzfile['arr_4']
Cr = npzfile['arr_5']


print(Y_dct[0][0])
print(Y[0][0])

print(Y.shape)
