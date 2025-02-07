import numpy as np
import scipy.sparse

from matplotlib import pyplot as plt
import csv
import sys

import scipy.interpolate as interpolate
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import UnivariateSpline

for case_path in np.arange(100, 600, 200):
    case_path = str(case_path) + '/'
    gain = np.load(case_path + 'gain_curve.npz')
    weight = []

    for S in range(len(gain["St"])):
        St = gain["St"][S]
        fq = np.load(case_path + 'fq_' + str("{:4f}".format(St)) + ".npz")
        temp = -gain["gain"][S] * gain["gain"][S] * np.imag(np.dot(fq["f"], fq["q"]))
        weight.append(temp)

    #weight = weight/ np.linalg.norm(weight)

    np.savez(case_path + 'spline_weight', St = gain["St"], weight = weight)

    weights = np.load(case_path + 'spline_weight.npz')
    gain = np.load(case_path + 'gain_curve.npz')

    St = gain["St"]
    gain_y = gain["gain"] 
    weight_values = weights["weight"]

    spline = UnivariateSpline(St, gain_y, w=weight_values, s=0)

    St_fine = np.linspace(St.min(), St.max(), 300)
    gain_fine = spline(St_fine)

    plt.plot(St_fine, gain_fine, color='blue', label='Weighted B-Spline Fit')
    plt.scatter(St, gain_y, color='red', label='Original Data', alpha=0.8)

#plt.xlim([0.0, 0.2])
#plt.figure(figsize=(8, 4))
plt.grid(True)
plt.yscale("log")
plt.ylabel("gain")
plt.xlabel("St")
plt.savefig("all_gain_w_spline.png")