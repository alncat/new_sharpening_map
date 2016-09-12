# coding: utf-8
import plotData
import pandas as pd
import numpy as np
import shutil
import os
new_sharp, sharp, pdbs = plotData.read_data('./testLog')
pList = pd.read_csv('pList', header=None)
#pList1 = pd.read_csv('pList1', header=None)
#pList = pd.concat([pList, pList1])
pdbs = np.array(pdbs)
diffs = np.setdiff1d(pList, pdbs)
#for diff in diffs:
#    if os.path.isdir(diff):
#        shutil.rmtree(diff)
import matplotlib.pyplot as plt
plt.figure()
fig, ax = plt.subplots()
import numpy as np
a_heights, a_bins = np.histogram(new_sharp)
b_heights, b_bins = np.histogram(sharp)
b_heights, b_bins = np.histogram(sharp, bins=a_bins)
width = (a_bins[1] - a_bins[0])/2
ax.bar(a_bins[:-1]+width, a_heights, width=width, facecolor='cornflowerblue', label='new sharp')
ax.bar(b_bins[:-1], b_heights, width=width, facecolor='seagreen', label='sharp')
plt.legend(loc=2)
plt.show()
