#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np

path="result2"

def ProcessTimeSize2AccselSize(path):
    data = np.loadtxt(path)

    num_workers = data[:, 0]
    time = data[:, 1]
    
    accsel = np.zeros(time.size)
    for i in range(time.size):
        accsel[i] = time[0] / time[i]

    return [num_workers, accsel]

size, time = ProcessTimeSize2AccselSize(path)

plt.scatter(size, time, color='blue', marker='x')

plt.grid()
plt.xlabel("Number threads")
plt.ylabel("Acceleration")

plt.show()
# plt.savefig(path + '.jpg')