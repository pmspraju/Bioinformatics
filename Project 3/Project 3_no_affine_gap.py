# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 22:15:15 2021

@author: pmspr
"""
# import packages
import os
import numpy as np

# Set the local paths for data
path = r'C:\Users\pmspr\Documents\HS\MS\Sem 4\EECS 730\Bioinformatics\Project 3\Docs'
f1 = os.path.join(path, 'file1.fasta')
f2 = os.path.join(path, 'file2.fasta')
seq = []
# Read file1.fasta
with open(f1) as fl1:
    for line in fl1:
        if(line[0].strip() != '>'):
            seq.append(line.strip())

# Read file2.fasta
with open(f2) as fl2:
    for line in fl2:
        if(line[0].strip() != '>'):
            seq.append(line.strip())

print('Sequence 1 - {}'.format(seq[0]))
print('Sequence 2 - {}'.format(seq[1]))

# Determine the lengths.
ls0 = len(seq[0]); ncols = ls0 + 1
ls1 = len(seq[1]); nrows = ls1 + 1

# Create numpy arrays for each input sequence
s0_arr = np.array(['*' for i in range(0, ncols)], dtype=np.object)
s0_arr[1:] = [seq[0][i] for i in range(0, ls0)]
#print(s0_arr)

s1_arr = np.array(['*' for i in range(0, nrows)], dtype=np.object)
s1_arr[1:] = [seq[1][i] for i in range(0, ls1)]
#print(s1_arr)

# Set the scoring parameters
match_score = 1.0
mismatch_score = -2.0
gap_score = -3

# Initialize the matrices with default values
sub_arr = np.zeros(shape=(nrows,ncols),dtype=np.int64)
sub_arr[0,1:] = np.linspace(gap_score, gap_score*ls0, ls0,dtype=np.int64)
sub_arr[1:,0] = np.linspace(gap_score, gap_score*ls1, ls1,dtype=np.int64)
trace_arr = np.array([['x' for i in range(0,ncols)]for j in range(0,nrows)], dtype=np.object)
trace_arr[0,1:] = np.array(['l' for i in range(0,ls0)],dtype=np.object)
trace_arr[1:,0] = np.array(['u' for i in range(0,ls1)],dtype=np.object)

# Traverse the matrices and using Needleman algorithm, fill the values
# substitution matrix(i,j) = max(i-1,j-1 & i,j-1 & i-1,j) 
# traceback matrix(i,j) = argmax(i-1,j-1 & i,j-1 & i-1,j)
neighbors = np.array([0,0,0], dtype=np.int64)
for i in range(1,nrows):
    for j in range(1,ncols):
        if(s1_arr[i] == s0_arr[j]):
            neighbors[0] = sub_arr[i-1,j-1] + match_score
        if(s1_arr[i] != s0_arr[j]):
            neighbors[0] = sub_arr[i-1,j-1] + mismatch_score
        neighbors[1] = sub_arr[i,j-1] + gap_score
        neighbors[2] = sub_arr[i-1,j] + gap_score
        sub_arr[i,j] = np.amax(neighbors)
        if(np.argmax(neighbors) == 0):
            trace_arr[i,j] = 'd'
        if(np.argmax(neighbors) == 1):
            trace_arr[i,j] = 'l'
        if(np.argmax(neighbors) == 2):
            trace_arr[i,j] = 'u'

print('Substitution Matrix:')
print(sub_arr)
print()
print('Trace back Matrix:')
print(trace_arr)

# Determine the alignment using traceback matrix
# Deter mine the score using scoring parameters
i = trace_arr.shape[0] - 1
j = trace_arr.shape[1] - 1
a0 = []; a1 = []; a2 = [];
score = 0
while(trace_arr[i,j] != 'x'):
    if(trace_arr[i,j] == 'd'):
        a0.append(s0_arr[j])
        a1.append(s1_arr[i])
        a2.append('|')
        score = score + match_score
        i = i-1;j=j-1;
    if(trace_arr[i,j] == 'l'):
        a0.append(s0_arr[j])
        a1.append('-')
        a2.append(' ')
        score = score + gap_score
        i=i;j=j-1;
    if(trace_arr[i,j] == 'u'):
        a0.append('-')
        a1.append(s1_arr[i])
        a2.append(' ')
        score = score + gap_score
        i=i-1;j=j;
    #print(i,j)
as0 = ''.join(a0[::-1]); as1 = ''.join(a1[::-1]); as2 = ''.join(a2[::-1])

# Pring the alignment to console
print()
print('Global Alignment without Affine gap with score = {}:'.format(score))
print(as0)
print(as2)
print(as1)