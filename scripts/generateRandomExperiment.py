import sys
import random

outDirectory = sys.argv[1] + '/'

rootLabel = '0'
numChildren = 4
numLevels = 3
queue = [rootLabel]
newNodes = []
numTrajectories = 324

for level in xrange(numLevels):
    print level
    numNodesInLevel = len(queue)
    clusters = [[] for i in xrange(numNodesInLevel)]
    for i in xrange(numTrajectories):
        clusterIndex = random.randint(0,numNodesInLevel - 1)
        clusters[clusterIndex].append(i)
    #write files
    print 'lengths',len(clusters),len(queue)
    count = 0
    for label in queue:
        f = open(outDirectory + label + '.txt','w')
        indicesInCluster = clusters[count]
        f.write('2\nheartrate\nspeed\n%d\n'%(len(indicesInCluster),))
        for i in indicesInCluster:
            f.write('%d\n'%(i,))
        count += 1
    #generate file names
    for node in queue:
        for cIndex in xrange(numChildren):
            cLabel = node + '_%d' % (cIndex,)
            newNodes.append(cLabel)
    
    #initialize
    queue = newNodes
    newNodes = []
