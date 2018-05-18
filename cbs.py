import sys
import numpy as np
import math
import da
from operator import itemgetter

"""
Implementation of task 2 in Computational Genomics assignment 2
Kinsey Reeves 18/05/2018
695705
"""

final_segs = []
Z_THRESH = 10

def setup(filename):
    
    f = open(filename, 'r')
    f.readline()
    data = []
    for i in f.readlines():
        line = i.strip('\n').strip('\r').split(' ')
        line = [line[0]] + [int(x) for x in line[1:]]
        data.append(line)

    #Get the first thirds median
    median = np.median([x[3] for x in data][0:int(len(data)*(1.0/3.0))])

    for i in range(0,len(data)):
        lr = (math.log(float(data[i][3])/float(median), 2))
        if lr>2 or lr<-5: lr = 0
        data[i][3] = lr
    return data

def cbs(data):
    sum_i = 0
    found = False
    best = -1*sys.maxsize
    best_i = 0
    best_j = 0

    n = len(data)
    sum_n = sum((x[3] for x in data))

    for i in range(0, n-1):
        sum_j = data[idx(i,n)][3]
        sum_j = 0

        for j in range(i+1, i+n):
            
            sum_j+=data[idx(j,n)][3]
            z = ( (1.0/(j-i)) + (1.0/(n-j+i)) )**(-1.0/2)
            z*= ( ((sum_j-sum_i)/(j-i)) - ((sum_n-sum_j+sum_i)/(n-j+i)))   
            z = abs(z)
        
            if(z>=Z_THRESH):
                if(z>best):
                    found = True
                    best = z
                    best_i = i
                    best_j = j
                
    if(found):
        segs = splice(best_i, best_j, data)
        final_segs.append(segs[0])
        cbs(segs[0])
        cbs(segs[1])
    else:
        return

#safe index
def idx(i, n):
    if(i >= n):
        return i-n
    else:
        return i

def splice(i, j, data):
    """
    safe splice
    input the full dataset and the 
    i and j coords and it will circularly 
    splice the segment and return two new 
    segments
    """
    n = len(data)
    out = ()
    splice_a = []
    splice_b = []
    i = idx(i, n)
    j = idx(j, n)
    
    if(i < j):
        splice_a = data[i+1:j+1]
        splice_b = data[0:i+1] + data[j+1:]
    else:
        splice_b = data[i+1:] + data[0:j+1]
        splice_a = data[j+1:i+1]

    return (splice_a, splice_b)
    
def print_segs(segs):

    if(segs):
        print("SEG A")
        for i in range(0,len(segs[0])):
            print(str(segs[0][i]) + " " + str(i) )
        print("SEG B")
        for i in range(0,len(segs[1])):
            print(str(segs[1][i]) + " " + str(i) )


data = setup("tumor.txt")
a = cbs(data)
out = []
for segment in final_segs:
    out.append(segment[0][0:2] + segment[-1][2:])

out = sorted(out, key=itemgetter(1))

for seg in out:
    print("{0}\t{1}\t{2}\t{3:.2f}".format(seg[0], seg[1], seg[2], seg[3]))
    





    




