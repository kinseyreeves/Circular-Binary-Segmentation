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
Z_THRESH = 8

def cut_contig_segs(seg, num):
    '''
    If a single segment isn't contiguous it cuts
    it into multiple. This is post-processing
    for the output
    '''
    out_segs = []
    if(len(seg)<=1):
        out_segs.append(seg)
        return output_segs
    step_size = seg[1][1]-seg[0][1]
    
    for i in seg:
        i.append(num)
    idx = 0
    for i in range(0, len(seg)-1):
        if(seg[i+1][1] - seg[i][1]!=step_size):
            out_segs.append(seg[idx:i+1])
            idx = i+1
    out_segs.append(seg[idx:])
    
    return out_segs

def setup(filename):
    '''
    Reads in the file and outputs a list of lists.
    Each item in the list is a row of the bin data.
    '''
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
    #Precomputes all sums 0 to n
    sums = []
    sum_i = 0
    
    for bi in data:
        sum_i+=bi[3]
        sums.append(sum_i)
    
    return sums

def cbs(data):
    '''
    1. Calculates the segment with highest Z score > thresh.
    If a segment is found it recursively does 1 again.
    Otherwise it outputs the segment to the global final_segs array.
    '''
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
    if(found):
        segs = splice(best_i, best_j, data)
        cbs(segs[0])
        cbs(segs[1])
    else:
        final_segs.append(data)


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

    return (splice_a, splice_b)
    
# Start of program

if(len(sys.argv) < 4):
    print("Usage: python cbs.py <input bins> <Z_thresh> <output file>")
    exit(0)
else:
    input_file = sys.argv[1]
    Z_THRESH = int(sys.argv[2])
    output_file = sys.argv[3]

output_file = open(output_file, "w")
data = setup("tumor.txt")
step_size = data[0][2] - data[0][1]

cbs(data)
output_segs = []

idx = 0
#We now need to cut segments so the indexes don't overlap
for j in final_segs:
    output_segs += cut_contig_segs(j, idx)
    idx+=1

output_segs.sort(key=lambda x: x[0][1])
out = []

#Get it into a format as given
for segment in output_segs:

    avg = sum([x[3] for x in segment])/len(segment)
    if(abs(avg)>0.1):
        out.append(segment[0][0:2] + [segment[-1][2]] + [avg] + [segment[0][-1]] + [len(segment)])

output_file.write("CHR\tSTART\tEND\tRATIO\tSEG NO.\n")
for seg in out:
    output_file.write("{0}\t{1}\t{2}\t{3:.2f}\t{4}\n".format(seg[0], seg[1], seg[2], seg[3], seg[4]))

output_file.close()



    




