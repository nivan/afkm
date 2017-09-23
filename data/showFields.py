'''
Show all different interpolation methods for imshow
'''

import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib import cm
import sys
from matplotlib.backends.backend_pdf import PdfPages

# from the docs:

def loadField(filename):
    reader = csv.reader(open(filename))
    field = []
    for line in reader:
        values = [float(t) for t in line]
        field.append(values)
    return field

def parseTrajData(filename):
    #
    f = open(filename)
    #read bbox
    bbox = [float(t)  for t in f.readline().strip().split(' ')]
    #
    line = f.readline().strip()
    numTrajectories = int(line.split(' ')[0])
    numAttributes   = int(line.split(' ')[1])

    trajectories = []
    for i in xrange(numTrajectories):
        numPoints = int(f.readline())
        traj = []
        for pIndex in xrange(numPoints):
            tokens = f.readline().strip().split(' ')
            point = [float(x) for x in tokens]
            traj.append(point)
        trajectories.append(traj)    
        
    #
    return (bbox, trajectories)

def parseClusterResultFile(filename):
    f = open(filename)
    result = []
    for line in f:
        result.append(int(line))
    return result


################################3333333

directory   = sys.argv[1]
trajFile    = 'synthetic.txt'
clusterFile = 'cluster_result_synthetic.txt'
groundTruthFile = 'generated_clusters.txt'

#
bbox, trajectories     = parseTrajData(directory + trajFile)
mapTrajectoryToCluster = parseClusterResultFile(directory + clusterFile)
groundTruthCluster     = parseClusterResultFile(directory + groundTruthFile)
numTrajectories = zip(trajectories,mapTrajectoryToCluster,groundTruthCluster)

#
numFields = 4
numAttributes = 4
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

for ax, field, title in zip(axes.flat, fields, titles):
    if title.find('Cluster') == -1:
        imgplot = ax.imshow(field, interpolation='none',cmap=cm.afmhot)
        ax.set_title(title)
        fig.colorbar(imgplot,ax = ax,shrink=0.9)
    else:
        ax.set_title(title)
        clusterIndex = int(title.split(' ')[1])
        x_coords = []
        y_coords = []
        for i in xrange(len(numTrajectories)):
            trajec  = numTrajectories[i]
            traj    = trajec[0]
            cluster = trajec[1]
            gtCluster = trajec[2]
            if cluster != clusterIndex:
                continue
            x_coords = [t[0] for t in traj]
            y_coords = [t[1] for t in traj]
            ax.plot(x_coords,y_coords,color=colors[cluster],linewidth=3)
            #ax.plot(x_coords,y_coords,color=colors[gtCluster],linewidth=3)

#
pp = PdfPages(directory + 'result_fig.pdf')
plt.savefig(pp, format='pdf')
pp.close()
#
fig.savefig(directory + 'foo.png',bbox_inches='tight')
plt.show()

