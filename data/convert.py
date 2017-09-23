import sys
import csv

filename = sys.argv[1]
reader = csv.reader(open(filename))

#
state = 0 #read experiment. 1 means read measurement

data = {}
numberOfSamples = -1
experiment_name = -1
count = 0
numMethods = 0
for line in reader:
    if state == 0:
        experiment_name = line[0]
        numberOfSamples = int(line[1])
        data[experiment_name] = []
        state = 1
        numMethods += 1
    elif state == 1:
        count += 1
        error = line[0]
        numIterations = line[1]
        elapsedTime = line[2]
        data[experiment_name].append((error,numIterations,elapsedTime))
        if count == numberOfSamples:
            state = 0
            count = 0

#print data
f = open(sys.argv[2],'w')
header = 'method,error,num_iter,time'
numExperiments = numberOfSamples
f.write(header + '\n')
methods = data.keys()

for method in data:
    methodData = data[method]
    for info in methodData:
        f.write('%s,%s,%s,%s\n'%(method,info[0],info[1],info[2]))
f.close()

            
# #print data
# f = open('out.csv','w')
# header = ''

# numExperiments = -1

# for method in data:
#     if header == '':
#         header += '%s_error,%s_num_iter,%s_time'%(method,method,method)
#         numExperiments = len(data[method])
#     else:
#         header += ',%s_error,%s_num_iter,%s_time'%(method,method,method)
#         assert(numExperiments == len(data[method]))
        
# f.write(header + '\n')
# methods = data.keys()
# for i in xrange(numExperiments):
#     line = ''
#     for mIndex in xrange(numMethods):
#         methodData = data[methods[mIndex]]
#         info = methodData[i]
#         line += '%s,%s,%s' % info
#         if mIndex < (numMethods-1):
#             line += ','
#     f.write(line + '\n')
# f.close()
