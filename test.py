# def idx(i, n):
#     if(i >= n):
#         return i-n
#     else:
#         return i

# print(2**(-(1/2)))

# print(1/(2**(1/2)))

#print(idx(22,10))
import math
import numpy as np


def setup(filename):
    
    f = open(filename, 'r')
    f.readline()
    data = []
    for i in f.readlines():
        line = i.strip('\n').strip('\r').split(' ')
        line = [line[0]] + [int(x) for x in line[1:]]
        data.append(line)

    median = np.median([x[3] for x in data][0:int(len(data)*(1.0/3.0))])

    for i in range(0,len(data)):
        lr = (math.log(float(data[i][3])/float(median), 2))
        if lr>2 or lr<-5: lr = 0
        data[i][3] = lr
    return data

a = setup("tumor.txt")

b= open("out.csv", "w")

for line in a:
    b.write(str(line[1]) + "," + str(line[3]) + "\n")
    