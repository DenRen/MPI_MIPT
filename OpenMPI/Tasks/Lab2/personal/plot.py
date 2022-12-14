#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np

path="result"

def ProcessTimeSize2AccselSize(path):
    data = np.loadtxt(path)

    num_workers = data[:, 0]
    time = data[:, 1]

    accsel = np.zeros(time.size)
    for i in range(time.size):
        accsel[i] = time[0] / time[i]

    return [num_workers, accsel]

def ProcessTimeSize2EffectiveSize(path):
    data = np.loadtxt(path)

    num_workers = data[:, 0]
    time = data[:, 1]

    eff = np.zeros(time.size)
    for i in range(time.size):
        accel = time[0] / time[i]
        eff[i] = accel / num_workers[i]

    return [num_workers, eff]

num_threads, acc = ProcessTimeSize2AccselSize(path)
num_threads, eff = ProcessTimeSize2EffectiveSize(path)

plt.subplot(2, 1, 1)
plt.plot(num_threads, acc, color='blue', marker='x')
plt.grid()
plt.xlabel("Number threads")
plt.ylabel("Acceleration")

plt.subplot(2, 1, 2)
plt.plot(num_threads, eff, color='blue', marker='x')
plt.grid()
plt.xlabel("Number threads")
plt.ylabel("Effective")

plt.show()
