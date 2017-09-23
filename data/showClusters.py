import sys
import random
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel

def lerp(v,w,factor):
    return (1-factor) * v + factor * w

def parseData(filename):
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

def parseResultFile(filename):
    f = open(filename)
    result = []
    for line in f:
        result.append(int(line))
    return result

if __name__ == '__main__':
    bbox, trajectories     = parseData(sys.argv[1])
    mapTrajectoryToCluster = parseResultFile(sys.argv[2])
    numTrajectories = zip(trajectories,mapTrajectoryToCluster)

    #
    x_coords = []
    y_coords = []
    #for trajec in numTrajectories:
    colors=['#e41a1c','#377eb8','#4daf4a','#984ea3']
    for i in xrange(len(numTrajectories)):
        trajec  = numTrajectories[i]
        traj    = trajec[0]
        cluster = trajec[1]
        if cluster != 0:
            continue
        #traj, cluster = numTrajectories[0]
        x_coords = [t[0] for t in traj]
        y_coords = [t[1] for t in traj]
        plt.plot(x_coords,y_coords,color=colors[cluster],linewidth=3)
    #plt.legend(colors , ('Cluster 1','Cluster 2','Cluster 3','Cluster 4'))
    # plt.tick_params(
    # axis='both',          # changes apply to the x-axis
    # which='both',      # both major and minor ticks are affected
    # bottom='off',      # ticks along the bottom edge are off
    # top='off',         # ticks along the top edge are off
    # left='off',
    # right='off',
    # labelbottom='off') # labels along the bottom edge are off
    plt.show()
