import sys

f = open(sys.argv[1])
outFile = open(sys.argv[2],'w')

#
gridLine = f.readline()
tokens = gridLine.split(' ')
tokens = [t.strip() for t in tokens]
outFile.write('%s %s %s %s %s %s\n' % (tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5]))
#
numAttribs = int(f.readline())
for i in xrange(numAttribs):
    f.readline()
#
numAttribs = int(f.readline())
for i in xrange(numAttribs):
    f.readline()
numAttribs -= 3
#
numTrajs = int(f.readline())
outFile.write('%d %d\n'%(numTrajs,numAttribs))
for i in xrange(numTrajs):
    #ignore statistics
    f.readline()
    lineNumPoints = f.readline()
    numPoints = int(lineNumPoints)
    outFile.write(lineNumPoints)
    for j in xrange(numPoints):
        outFile.write(f.readline())
        
