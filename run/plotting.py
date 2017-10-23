import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams


image = np.load('../data/outfiles/1996-05-16_18_OBJECT_1.0.npy')
land  = np.load('../data/land.npy')

v    = ma.masked_equal(image[1], 0)
u    = ma.masked_equal(image[2], 0)
ivt  = ma.masked_equal(image[0], 0)
path = ma.masked_equal(image[3], 0)




fig, ax = plt.subplots(1,1)
ax.imshow(land, alpha = .5)
ax.imshow(ivt, alpha=.5)
ax.imshow(path, alpha=.5)

#plt.quiver(v,u, scale=1000., alpha = .5)
plt.show()