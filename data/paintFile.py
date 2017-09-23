'''
Show all different interpolation methods for imshow
'''

import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib import cm
import sys
from matplotlib.backends.backend_pdf import PdfPages

#
numFields = 1
numAttributes = 1
fields = []
titles = []
for fIndex in xrange(numFields):
    for aIndex in xrange(numAttributes):
        if aIndex < (numAttributes - 1):
            fields.append(loadField('%s/field_%d_%d.txt'%(directory,fIndex,aIndex)))
            titles.append('Field %d Att %d' % (fIndex,aIndex))
        else:
            fields.append(('cluster',fIndex))
            titles.append('Cluster %d' % (fIndex,))

fig, axes = plt.subplots(numFields, numAttributes, figsize=(12, 6),
                         subplot_kw={'xticks': [], 'yticks': []})
fig.subplots_adjust(hspace=0.3, wspace=0.05)
colors=['#e41a1c','#377eb8','#4daf4a','#984ea3']

for axe in axes.flat:
    imgplot = ax.imshow(field, interpolation='none',cmap=cm.afmhot)
    ax.set_title('MyPlot')
    fig.colorbar(imgplot,ax = ax,shrink=0.9)

