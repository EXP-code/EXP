import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt

if not os.path.exists(sys.argv[1]):
    print("File <{}> does not exist?".format(sys.argv[1]))
    exit(1)

file = open(sys.argv[1], 'rb')
db = pickle.load(file)
file.close()

image = db['image']
lower = db['lower']
upper = db['upper']
ngrid = db['ngrid']

for v in image:
    im = plt.imshow(image[v]["xy"].transpose(),
                    extent=[lower[0], upper[0], lower[1], upper[1]])
    plt.colorbar(im)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Time={}'.format(v))
    plt.show()
