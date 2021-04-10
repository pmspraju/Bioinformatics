# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 22:11:49 2021

@author: pmspr
"""
# import packages
import os
import numpy as np

# Set the local paths for data
path = r'C:\Users\pmspr\Documents\HS\MS\Sem 4\EECS 730\Bioinformatics\Project 3\Docs'
f1 = os.path.join(path, 'file1.fasta') # Change the file name accordingly.
f2 = os.path.join(path, 'file2.fasta') # Change the file name accordingly.
seq = []
#Read file1.fasta
with open(f1) as fl1:
    for line in fl1:
        if(line[0].strip() != '>'):
            seq.append(line.strip())

#Read file2.fasta
with open(f2) as fl2:
    for line in fl2:
        if(line[0].strip() != '>'):
            seq.append(line.strip())

print('Sequence 1 - {}'.format(seq[0]))
print('Sequence 2 - {}'.format(seq[1]))

# Determine the lengths
ls0 = len(seq[0]); ncols = ls0 + 1
ls1 = len(seq[1]); nrows = ls1 + 1

# Create numpy arrays from the strings
s0_arr = np.array(['*' for i in range(0, ncols)], dtype=np.object)
s0_arr[1:] = [seq[0][i] for i in range(0, ls0)]
#print(s0_arr)

s1_arr = np.array(['*' for i in range(0, nrows)], dtype=np.object)
s1_arr[1:] = [seq[1][i] for i in range(0, ls1)]
#print(s1_arr)

# Define the scoring parameters
match_score = 1.0
mismatch_score = -2.0
gap_open = -5
gap_extn = -1

# Define the upper layer matrix and initialize
# Initialize first column = gap_open + row_index * gap_extn
upp_arr = np.zeros(shape=(nrows,ncols),dtype=np.int64)
upp_arr[0,0] = -1000
upp_arr[1:,0] = np.array([gap_open + i*gap_extn for i in range(1,nrows)],dtype=np.int64)
#print(upp_arr)

# Define the left layer matrix and initialize
# Initialize first row = gap_open + col_index * gap_extn
lft_arr = np.zeros(shape=(nrows,ncols),dtype=np.int64)
lft_arr[0,0] = -1000
lft_arr[0,1:] = np.array([gap_open + i*gap_extn for i in range(1,ncols)],dtype=np.int64)
#print(lft_arr)

# Define the diagonal layer matrix and initialize
# Initialize first column, row = 0
dig_arr = np.zeros(shape=(nrows,ncols),dtype=np.int64)
dig_arr[0,0] = -1000
#print(dig_arr)

# Define the diagonal layer matrix and initialize
# Initialize first column = gap_open + row_index * gap_extn
# Initialize first row = gap_open + col_index * gap_extn
sub_arr = np.zeros(shape=(nrows,ncols),dtype=np.int64)
sub_arr[0,1:] = np.array([gap_open + i*gap_extn for i in range(1,ncols)],dtype=np.int64)
sub_arr[1:,0] = np.array([gap_open + i*gap_extn for i in range(1,nrows)],dtype=np.int64)
#print(sub_arr)


# Define trace back matrix
# Initialize first column = 'u'
# Initialize first row = 'l'
trace_arr = np.array([['x' for i in range(0,ncols)]for j in range(0,nrows)], dtype=np.object)
trace_arr[0,1:] = np.array(['l' for i in range(0,ls0)],dtype=np.object)
trace_arr[1:,0] = np.array(['u' for i in range(0,ls1)],dtype=np.object)

# Fill the intermediate arrays
neighbors = np.array([0,0], dtype=np.int64)
fneigh = np.array([0,0,0], dtype=np.int64)
for i in range(1,nrows):
    for j in range(1,ncols):
        
        # Fill upper matrix
        neighbors = [0,0]
        neighbors[0] = sub_arr[i,j-1] + gap_open + gap_extn
        neighbors[1] = upp_arr[i,j-1] + gap_extn
        upp_arr[i,j] = np.amax(neighbors)
        
        # Fill left matrix
        neighbors = [0,0]
        neighbors[0] = sub_arr[i-1,j] + gap_open + gap_extn
        neighbors[1] = lft_arr[i-1,j] + gap_extn
        lft_arr[i,j] = np.amax(neighbors)
        
        # fill diagonal matrix
        if(s1_arr[i] == s0_arr[j]):
            dig_arr[i,j] = sub_arr[i-1,j-1] + match_score
        if(s1_arr[i] != s0_arr[j]):
            dig_arr[i,j] = sub_arr[i-1,j-1] + mismatch_score
        
        # fill traceback matrix
        fneigh = [0,0,0]
        fneigh = np.array([dig_arr[i,j],upp_arr[i,j],lft_arr[i,j]],dtype=np.int64)
        sub_arr[i,j] = np.amax(fneigh)
        if(np.argmax(fneigh) == 0):
            trace_arr[i,j] = 'd'
        if(np.argmax(fneigh) == 1):
            trace_arr[i,j] = 'l'
        if(np.argmax(fneigh) == 2):
            trace_arr[i,j] = 'u'

print('Upper layer Matrix:')
print(upp_arr)
print()

print('Left layer Matrix:')
print(lft_arr)
print()

print('Diagonal layer Matrix:')
print(dig_arr)
print()

print('Substitution Matrix:')
print(sub_arr)
print()

print('Trace back Matrix:')
print(trace_arr)

# Determine the alignment using traceback matrix
i = trace_arr.shape[0] - 1
j = trace_arr.shape[1] - 1
a0 = []; a1 = []; a2 = [];
score = 0;ugap_ind = 0; lgap_ind = 0;

# Traverse from right bottom to the top left to create the alignment
while(trace_arr[i,j] != 'x'):
    if(trace_arr[i,j] == 'd'):
        a0.append(s0_arr[j])
        a1.append(s1_arr[i])
        if(s0_arr[j] == s1_arr[i]):
            a2.append('|')
            score = score + match_score
            
        if(s0_arr[j] != s1_arr[i]):
            a2.append('*')
            score = score + mismatch_score
        i = i-1;j=j-1;
    if(trace_arr[i,j] == 'l'):
        a0.append(s0_arr[j])
        a1.append('-')
        a2.append(' ')
        lgap_ind = lgap_ind + 1
        score = score + gap_open + lgap_ind*gap_extn
        i=i;j=j-1;
    if(trace_arr[i,j] == 'u'):
        a0.append('-')
        a1.append(s1_arr[i])
        a2.append(' ')
        ugap_ind = ugap_ind + 1
        score = score + gap_open + ugap_ind*gap_extn
        i=i-1;j=j;
    
as0 = ''.join(a0[::-1]); as1 = ''.join(a1[::-1]); as2 = ''.join(a2[::-1])

# print the output to the console
print()
print('Global Alignment without Affine gap with score = {}:'.format(score))
print(as0)
print(as2)
print(as1)