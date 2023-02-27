# find the distribution of combo from a demo image

import numpy as np
from osgeo import gdal

img = gdal.Open(r"result\best_combo\ca_mesma_demo.dat")

im_width = img.RasterXSize
im_height = img.RasterYSize
im_band = img.RasterCount
im_data = img.ReadAsArray(0, 0, im_width, im_height)
first = im_data[0, :, :]
number, counts = np.unique(first, return_counts=True)
output = np.zeros((len(number), 3))
output[:, 0] = number.astype(int)
output[:, 1] = counts
output[:, 2] = output[:, 1] / np.sum(output[:, 1])

out = output[np.argsort(-output[:, 1])]

np.savetxt('result\best_combo\mesma_stats.txt', out, fmt='%d  %d  %5.5f')
