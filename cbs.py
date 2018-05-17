import sys
import numpy as np
import math

"""
Implementation of task 2 in Computational Genomics assignment 2
"""


Z_THRESH = 1

def setup():
    
    f = open("sample.txt")
    f.readline()
    data = []
    for i in f.readlines():
        line = i.strip('\n').strip('\r').split(' ')
        line = [line[0]] + [int(x) for x in line[1:]]
        data.append(line)

    #Get the first thirds median
    median = np.median([x[3] for x in data][0:int(len(data)*(1.0/3.0))])

    for i in range(0,len(data)):
        lr = math.log(float(data[i][3])/float(median), 2)
        if lr>2 or lr<-5: lr = 0
        data[i][3] = lr
    return data
    

def cbs(data, i_start, j_start, n, sum_n):
    sum_i = 0
    best = -99999
    best_i = 0
    best_j = 0
    for i in range(i_start, n):
        
        sum_j = 0
        sum_i+=data[i][3]

        for j in range(i+1, i+n):
            print(str(i) + ' full: ' + str(j) + ' actual: ' + str(idx(j,n)))

            sum_j+=data[idx(j,n)][3]
            z = ((1/(j-i)) + (1/(n-j+i)))**(-1/2)
            z*= ((sum_j-sum_i)/(j-i))-((sum_n-sum_j+sum_i)/(n-j+i))

            if(z>Z_THRESH):
                if(z>best):
                    best = z
                    best_i = i
                    best_j = j
                
    print(z)
    print(best_i)
    print(best_j)
    return splice(best_i, best_j, data)

#safe index
def idx(i, n):
    n_len = n
    if(i >= n_len):
        return i-n_len
    else:
        return i

#safe splice
#input the full dataset and the 
#i and j coords 
def splice(i, j, data):
    i+=1
    j+=1
    n = len(data)
    out = ()
    splice_a = []
    splice_b = []
    if(j < n):
        splice_a = data[i:j]
        splice_b = data[j:] + data[0:i]
        out = (splice_a, splice_b)
    else:
        splice_a = data[i:] + data[0:idx(j, n)]
        splice_b = data[idx(j, n):i]
        out = (splice_a, splice_b)
    return out
    

    




data = setup()
n = len(data)
sum_n = sum([x[3] for x in data])
# print(sum_n)


# print(n)
# print()

# print(idx(51, n))

    
# a = splice(40,55,data)[0]
# b = splice(40,55,data)[1]

a = cbs(data, 0, 1, len(data), sum_n)

for i in a[0]:
    print(i)
print()
for i in a[1]:
    print(i)





    




