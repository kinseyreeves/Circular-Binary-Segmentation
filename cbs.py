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
Z_THRESH = 5

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


def get_starts_ends(segment):
    return 

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
        
            if(z>=Z_THRESH and z>best):
                found = True
                best = z
                best_i = i
                best_j = j
                
    if(found):
        segs = splice(best_i, best_j, data)
        #print("here")
        #print_segs(segs)
        #a = input()
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
    out = ()
    splice_a = []
    splice_b = []
    i = idx(i, n)
    j = idx(j, n)
    i+=1
    j+=1
    if(i < j):
        splice_a = data[i:j]
        splice_b = data[j:] + data[0:i]
    else:
        splice_a = data[i:] + data[0:j]
        splice_b = data[j:i]

    return (splice_a, splice_b)
    

data = setup("tumor.txt")
a = cbs(data)
out = []
for segment in final_segs:
    avg = sum([x[3] for x in segment])/len(segment)
    #if(abs(avg)>0.1):
    out.append(segment[0][0:2] + [segment[-1][2]] + [avg] + [len(segment)])
    
    for line in segment:
        print(line)
    print("\n\n\n\n\n"*30)


out = sorted(out, key=itemgetter(1))

for seg in out:
    
    print("{0}\t{1}\t{2}\t{3:.2f}\t+{4}".format(seg[0], seg[1], seg[2], seg[3], seg[4]))
    





    




