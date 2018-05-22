import sys
import numpy as np
import math
import da
from operator import itemgetter


"""
Implementation of task 2 in Computational Genomics assignment 2
Kinsey Reeves 18/05/2018. 
695705
"""


final_segs = []
Z_THRESH = 1

def print_segs(segs):
    
    if(segs):
        print("SEG A")
        for i in range(0,len(segs[0])):
            print(str(segs[0][i]) + " " + str(i) )
        print("SEG B")
        for i in range(0,len(segs[1])):
            print(str(segs[1][i]) + " " + str(i) )

def print_seg(segs):
    for i in segs:
        print(i)


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


def precompute_sums(data):
    sums = []
    sum_i = 0
    
    for bi in data:
        #print(bi[3])
        sum_i+=bi[3]
        sums.append(sum_i)
    
    #print(sums)
    return sums

def cbs(data):
    sum_i = 0
    found = False
    best = -1*sys.maxsize
    best_i = 0
    best_j = 0
    sums = precompute_sums(data)

    n = len(data)
    sum_n = sum((x[3] for x in data))

    for i in range(0, n-1):
        sum_i = sums[i]
        sum_j = 0

        for j in range(i+1, n):
            
            sum_j = sums[j]
            z = ( (1.0/(j-i)) + (1.0/(n-j+i)) )**(-0.5)
            z*= ( ((sum_j-sum_i)/(j-i)) - ((sum_n-sum_j+sum_i)/(n-j+i)))   
            z = abs(z)
        
            if(z>=Z_THRESH):

                found = True
                if(z>best):
                    best = z
                    best_i = i
                    best_j = j

    #exit(0)            
    if(found):
        segs = splice(best_i, best_j, data)
        cbs(segs[0])
        cbs(segs[1])
    else:
        
        final_segs.append(data)
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
    splice_a = []
    splice_b = []
    i+=1
    j+=1
    splice_a = data[i:j]
    splice_b = data[j:] + data[0:i]
    print(data)
    exit(0)
    

    return (splice_a, splice_b)
    

data = setup("sample.txt")

a = cbs(data)


out = []
for segment in final_segs:

    avg = sum([x[3] for x in segment])/len(segment)
    if(abs(avg)>0.1):
        out.append(segment[0][0:2] + [segment[-1][2]] + [avg] + [len(segment)])
    
    for line in segment:
        print(line)
    print("\n\n\n\n\n"*30)


out = sorted(out, key=itemgetter(1))

for seg in out:
    
    print("{0}\t{1}\t{2}\t{3:.2f}\t+{4}".format(seg[0], seg[1], seg[2], seg[3], seg[4]))
    





    




