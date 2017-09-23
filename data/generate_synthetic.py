import sys
import random
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel
import math

#######################
# Auxiliary Functions #
#######################

def toDegree(x):
    return 180.0 * x / math.pi

def toRadians(x):
    return math.pi * x / 180.0

def directionToAngle(x,y):
    if x == 0:
        if y > 0:
            return 90
        else:
            return 270
    else:
        return toDegree(math.atan(y / x))


def lerp(v,w,factor):
    return (1-factor) * v + factor * w

def parseData(filename):
    #
    f = open(filename)
    #read bbox
    bbox = f.readline().strip()
    #read trajectories
    trajectories = []
    #
    currentTrajectory = []
    for line in f:
        if line == '0.0 0.0 0.0\n':
            trajectories.append(currentTrajectory)
            currentTrajectory = []
        else:
            currentTrajectory.append(line.strip())

    #
    return (bbox, trajectories)

########################
# Convertion Functions #
########################
            
def convertFileWebFormat(newFileName,bbox,trajectories):
    #v1
    fields  = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    #v2
    #fields  = [[-0.5,1,1,-0.5],[1,0.5,0.5,1],[1,-.3,-.3,1],[-.2,.5,.5,-.2]]
    #v3
    #fields  = [[-0.5,.8,-.8,-0.5],[1,-0.5,1,-1],[-.3,1,-.3,1],[.2,-.5,.5,.2]]
    #non overlapping fields
    #fields  = [[-1,-.8,-.8,-1],[-.5,0,0,-.5],[.3,.5,.5,.3],[.7,1,1,.7]]
    
    #
    xmin,xmax,ymin,ymax,tmin,tmax = bbox.split(' ')
    numSpatialBBox = [float(xmin),float(xmax),float(ymin),float(ymax)]

    #
    config = {'write_tgx':True, 'write_tgy':True, 'write_field_value':True}
    numAttributes = 0
    for key in config:
        if config[key]:
            numAttributes += 1
            
    #build header
    extraHeader = ''
    if config['write_tgx']:
        extraHeader += (' -1 1')
    if config['write_tgy']:
        extraHeader += (' -1 1')
    if config['write_field_value']:
        extraHeader += (' -1 1')
    
    #
    parameterConfig = '4\nlatitude,Latitute,number\nlongitude,Longitude,number\ntime,Time,number\ncluster,Cluster,number\n%d\nlatitude,Latitude,timeserie\nlongitude,Longitude,timeserie\ntime,Time,timeserie\n' % (numAttributes+3,)
    if config['write_tgx']:
        parameterConfig += ('tgX,tgX,timeserie\n')
    if config['write_tgy']:
        parameterConfig += ('tgY,tgY,timeserie\n')
    if config['write_field_value']:
        parameterConfig += ('field_value,field_value,timeserie\n')    

    #
    numberTrajectories = len(trajectories)

    #
    f = open(newFileName,'w')
    f.write(bbox + extraHeader + '\n' + parameterConfig + ('%d\n'%(numberTrajectories,)))
    
    #
    countsPerCluster = {}
    for i in xrange(len(fields)):
        countsPerCluster[i] = 0

    numTrajectories = []
    for i in xrange(len(trajectories)):
        if i % 2 == 0:
            #top center
            fieldIndex = random.randint(0,1)
        else:
            #bottom center
            fieldIndex = random.randint(2,3)        

        #
        f.write('0 0 0 %d\n' % fieldIndex)
        #
        traj = trajectories[i]
        numPoints = len(traj)
        f.write('%d\n'%(numPoints,))
        #
        countsPerCluster[fieldIndex] += 1
        numTraj = []
        #for point in traj:
        numPoints = len(traj)
        for j in xrange(numPoints):
            point = traj[j]
            x,y,t = point.split(' ')
            #

            if numPoints == 1:
                tgx = 0
                tgy = 0
            else:                
                numX = float(x)
                numY = float(y)
                numT = float(t)
                if j == (numPoints - 1):
                    otherPointX,otherPointY,otherPointT = [float(z) for z in traj[j-1].split(' ')]
                    tgx = (numX - otherPointX) / (numT - otherPointT)
                    tgy = (numY - otherPointY) / (numT - otherPointT)
                else:
                    otherPointX,otherPointY,otherPointT = [float(z) for z in traj[j+1].split(' ')]
                    tgx = (otherPointX - numX) / (otherPointT - numT)
                    tgy = (otherPointY - numY) / (otherPointT - numT)
            #
            value = getSpatialFactor(numSpatialBBox,(float(x),float(y)),fields[fieldIndex])
            cluster = (1.0*fieldIndex)/len(fields)

            #
            extraCoords = ''
            if config['write_tgx']:
                extraCoords += (' %f' % (tgx,))
            if config['write_tgy']:
                extraCoords += (' %f' % (tgy,))
            if config['write_field_value']:
                extraCoords += (' %f' % (value,))

            extraCoords += '\n'            
            f.write(point + extraCoords)            
            numTraj.append((x,y))
        #
        numTrajectories.append((numTraj,fieldIndex))

    #
    f.close()
    print countsPerCluster
    return numTrajectories

def getSpatialFactor(numBBox, point, field):
    #
    xmin,xmax,ymin,ymax = numBBox    
    x,y = point
    #
    lambdaX = (x - xmin)/(xmax - xmin)
    lambdaY = (y - ymin)/(ymax - ymin)

    #
    v0,v1,v2,v3 = field
    #
    return lerp(lerp(v0,v1,lambdaX),lerp(v3,v2,lambdaX),lambdaY)

def radialField(numBBox, point):
    #
    xmin,xmax,ymin,ymax = numBBox    
    x,y = point
    #
    lambdaX = (x - xmin)/(xmax - xmin)
    lambdaY = (y - ymin)/(ymax - ymin)
    #
    return 2*( (lambdaX - 0.5) ** 2 + (lambdaY - 0.5) ** 2 )

def reverseRadial(numBBox, point):
    x = radialField(numBBox, point)
    return 1 - x

def getFieldEvaluationFunction(field):
    return lambda x,y:getSpatialFactor(x,y,field)

def mapTangent(numBBox, point):
    return directionToAngle(point[0],point[1]) / 360.0
        

def convertFile(newFileName,directory,bbox,trajectories):
    #v1
    fields  = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]    
    #v2
    #fields  = [[1,1,0,0],[1,0,0,1],[1,0,1,0],[0,1,0,1]]
    #v3
    #fields  = [[-0.5,.8,-.8,-0.5],[1,-0.5,1,-1],[-.3,1,-.3,1],[.2,-.5,.5,.2]]
    #non overlapping fields
    #fields  = [[-1,-.8,-.8,-1],[-.5,0,0,-.5],[.3,.5,.5,.3],[.7,1,1,.7]]    
    fieldFunctions = []
    for field in fields:
        fieldFunction = getFieldEvaluationFunction(field)
        fieldFunctions.append(fieldFunction)
 

    fieldFunctions = [getFieldEvaluationFunction([1,0,0,0]),radialField,getFieldEvaluationFunction([0,0,1,0]),reverseRadial]

    #
    fieldFunctions = [getFieldEvaluationFunction([1,0,0,0]),mapTangent,lambda x,y: mapTangent(x,[-t for t in y]),getFieldEvaluationFunction([0,0,1,0])]

    #
    xmin,xmax,ymin,ymax,tmin,tmax = bbox.split(' ')
    numSpatialBBox = [float(xmin),float(xmax),float(ymin),float(ymax)]

    #
    clusterFile = open(directory + '/generated_clusters.txt','w')
    #
    f = open(newFileName,'w')
    f.write(bbox + '\n')
    numAttributes = 3
    f.write('%d %d\n'%(len(trajectories),numAttributes)) # fixed one attribute for now
    
    #
    countsPerCluster = {}
    for i in xrange(len(fields)):
        countsPerCluster[i] = 0

    numTrajectories = []
    countTrajs = len(trajectories)
    for i in xrange(countTrajs):
        if i < (countTrajs/2):
            fieldIndex = (i%2)
        else:
            fieldIndex = (i%2) + 2
        # if i % 2 == 0:
        #     #top center
        #     fieldIndex = random.randint(0,1)
        # else:
        #     #bottom center
        #     fieldIndex = random.randint(2,3)
        traj = trajectories[i]
        numPoints = len(traj)
        clusterFile.write('%d\n' % (fieldIndex,))
        f.write('%d\n'%(numPoints,))
        #
        countsPerCluster[fieldIndex] += 1
        numTraj = []
        #for point in traj:
        numPoints = len(traj)
        for j in xrange(numPoints):
            point = traj[j]
            x,y,t = point.split(' ')
            #
            numX = float(x)
            numY = float(y)
            numT = float(t)

            if numPoints == 1:
                tgx = 0
                tgy = 0
            else:                
                if j == (numPoints - 1):
                    otherPointX,otherPointY,otherPointT = [float(z) for z in traj[j-1].split(' ')]
                    tgx = (numX - otherPointX) / (numT - otherPointT)
                    tgy = (numY - otherPointY) / (numT - otherPointT)
                else:
                    otherPointX,otherPointY,otherPointT = [float(z) for z in traj[j+1].split(' ')]
                    tgx = (otherPointX - numX) / (otherPointT - numT)
                    tgy = (otherPointY - numY) / (otherPointT - numT)
                                
            #
            #value = getSpatialFactor(numSpatialBBox,(float(x),float(y)),fields[fieldIndex])
            if fieldIndex == 1 or fieldIndex == 2:
                value = fieldFunctions[fieldIndex](numSpatialBBox,(tgx,tgy))
            else:                
                value = fieldFunctions[fieldIndex](numSpatialBBox,(numX,numY))

            cluster = (1.0*fieldIndex)/len(fields)
            if numAttributes == 2:
                #f.write(point + ' %f %f\n'%(value,cluster))
                f.write(point + ' %f %f\n'%(tgx,tgy))            
            elif numAttributes == 1:
                f.write(point + ' %f\n'%(value,))            
            elif numAttributes == 3:
                f.write(point + ' %f %f %f\n'%(value,tgx,tgy))            
            numTraj.append((x,y))
        #
        numTrajectories.append((numTraj,fieldIndex))

    #
    f.close()
    print countsPerCluster
    countsFile = open(outDirectory + '/generatedCounts.txt','w')
    countsFile.write(str(countsPerCluster))
    countsFile.close()
    return numTrajectories

if __name__ == '__main__':
    
    trajFile     = sys.argv[1]
    outDirectory = sys.argv[2]
    outFileName  = sys.argv[3]
    bbox, trajectories = parseData(trajFile)
    #
    # print 'bbox',bbox
    # print 'num trajectories %d' % (len(trajectories),)
    #
    #numTrajectories = convertFile('experiment_synthetic_weights/synthetic.txt','experiment_synthetic_weights/',bbox,trajectories)
    numTrajectories = convertFile(outDirectory + '/' + outFileName,outDirectory,bbox,trajectories)
    #numTrajectories = convertFileWebFormat('test_file.txt',bbox,trajectories)

    # #
    # x_coords = []
    # y_coords = []
    # #for trajec in numTrajectories:
    # colors=['#e41a1c','#377eb8','#4daf4a','#984ea3']
    # for i in xrange(len(numTrajectories)):
    #     trajec  = numTrajectories[i]
    #     traj    = trajec[0]
    #     cluster = trajec[1]
    #     #traj, cluster = numTrajectories[0]
    #     x_coords = [t[0] for t in traj]
    #     y_coords = [t[1] for t in traj]
    #     plt.plot(x_coords,y_coords,color=colors[cluster],linewidth=3)
    # plt.show()
