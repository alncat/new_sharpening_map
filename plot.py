# coding: utf-8
import plotData
import pandas as pd
import numpy as np
import shutil
import os
new_sharp, sharp, pdbs = plotData.read_data('./testLog')
pList = pd.read_csv('pList', header=None)
pList1 = pd.read_csv('pList1', header=None)
pList2 = pd.read_csv('pList2', header=None)
pList = pd.concat([pList, pList1, pList2])
stats = pd.read_csv('stat_log', header=None, sep=' ')
high_b = stats[2] > 4.0
pdbs = np.array(pdbs)
intsect = np.in1d(pdbs, stats[0][high_b])
new_sharp_int = []
sharp_int = []
print intsect.shape, len(pdbs)
for i in xrange(len(intsect)):
  if(intsect[i]):
    new_sharp_int.append(new_sharp[i])
    sharp_int.append(sharp[i])
print len(new_sharp_int)
diffs = np.setdiff1d(pList, pdbs)
#for diff in diffs:
#    if os.path.isdir(diff):
#        shutil.rmtree(diff)
import matplotlib.pyplot as plt
plt.figure()
fig, ax = plt.subplots()
import numpy as np
a_heights, a_bins = np.histogram(new_sharp_int)
#b_heights, b_bins = np.histogram(sharp_int)
b_heights, b_bins = np.histogram(sharp_int, bins=a_bins)
width = (a_bins[1] - a_bins[0])/2
ax.bar(a_bins[:-1]+width, a_heights, width=width, facecolor='cornflowerblue', label='new sharp')
ax.bar(b_bins[:-1], b_heights, width=width, facecolor='seagreen', label='sharp')
plt.legend(loc=2)
plt.show()
