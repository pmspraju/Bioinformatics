# -*- coding: utf-8 -*-
"""
Created on Thu May 13 16:58:06 2021

@author: pmspr
"""

# import packages
import os
import numpy as np
import pandas as pd
import math
import sys

# Print versions
np.set_printoptions(threshold=sys.maxsize)

# Take file names as input. # Change the file name accordingly.
ffasta = 'AAA.fa' 
fhmm = 'AAA.hmm'

# Set the local paths for data
path = r'C:\Users\pmspr\Documents\HS\MS\Sem 4\EECS 730\Bioinformatics\Project 4\Docs'
fasta = os.path.join(path, ffasta) 
hmm   = os.path.join(path, fhmm)
test  = os.path.join(path,'text.csv')

#Read hmm file in to dataframes
states = [];tran_cols = [];
match = []; insert = []; trans = [];
hmmoffset = 0; holdstate = ' ';
state_cnt = 0;lineoffset = 0;

with open(hmm) as flhmm:
    for line in flhmm:
        if (line.strip()[0:3] == '//'):
            break
            
        if (line[0:4].strip() == 'HMM'):
            states = [s for s in line.split(' ') if (s not in ['','HMM','\n'])]
            state_cnt = len(states)
            hmmoffset += 1
            continue 
        if (hmmoffset == 1):
            tran_cols = [s for s in line.strip('\n').split(' ') if (s not in ['','\n'])]
            hmmoffset += 1
            continue
        if (hmmoffset == 2):
            lineoffset += 1
            if (lineoffset == 1):
                match_row = []
                match_row = [s for s in line.strip('\n').split(' ') if (s not in ['','\n'])][0:state_cnt+1]
                holdstate = match_row[0]
                match.append(match_row)
            if (lineoffset == 2):
                insert.append([holdstate]+[s for s in line.strip('\n').split(' ') if (s not in ['','\n'])][0:state_cnt+1])
            if (lineoffset == 3):
                trans.append([holdstate]+[s for s in line.strip('\n').split(' ') if (s not in ['','\n'])][0:len(tran_cols)+1])
                lineoffset = 0; holdstate = ' ';
            
match_df = pd.DataFrame(match, columns=['state']+states); match_df['state'].iloc[0] = 0;
insert_df = pd.DataFrame(insert, columns=['state']+states); insert_df['state'].iloc[0] = 0;
trans_df =  pd.DataFrame(trans,columns=['state']+tran_cols); trans_df['state'].iloc[0] = 0;
trans_df.loc[trans_df['state'] == 0] = 1;

# Display sample dataframes
print(match_df.head())
print('Match emission probabilities')
print(insert_df.tail())
print('Insert emission probabilities')
print(trans_df.head())
print('Transition probabilities')

# Read the sequence file
seql = []
with open(fasta) as flfst:
    for line in flfst:
        if (line.strip()[0] != '>'):
            seq = line.strip('\n')
seql = [seq[s] for s in range(len(seq))]
print('Input sequence with length {}:\n{}'.format(len(seq),seq))

# Calculate matrices
ncols = match_df.shape[0] - 1     #Ignore compo
nrows = len(seq) 

# Declare and initialize

# Match matrix
vm = np.zeros(shape=(nrows+1,ncols+1),dtype=np.float64)
vm[0,0] = 0; vm[0,1:] = -1000; vm[1:,0] = -1000
vmt = np.array([['x' for i in range(0,ncols+1)]for j in range(0,nrows+1)], dtype=np.object)

# Insert matrix
vi = np.zeros(shape=(nrows+1,ncols+1),dtype=np.float64)
vit = np.array([['x' for i in range(0,ncols+1)]for j in range(0,nrows+1)], dtype=np.object)
vi[0:,0] = -1000

# Delete matrix
vd = np.array([0,0], dtype=np.float64)
vd = np.zeros(shape=(nrows+1,ncols+1),dtype=np.float64)
vdt = np.array([['x' for i in range(0,ncols+1)]for j in range(0,nrows+1)], dtype=np.object)
vd[0,0:] = -1000

# Terminate matrix
vt = np.zeros(shape=(nrows+1,ncols+1),dtype=np.float64)

# Trace matrix
tb = np.array([['x' for i in range(0,ncols+1)]for j in range(0,nrows+1)], dtype=np.object)

# Fill the matrice
n = 0.0; n1 = 0.0; n2 = 0.0; n3 = 0.0;
for j in range(1, nrows+1):
    for i in range(1,ncols+1):
        
        # Match matrix
        if (trans_df['m->m'].iloc[j-1] == '*'):
            trans_df['m->m'].iloc[j-1] = -1000
        
        if (trans_df['i->m'].iloc[j-1] == '*'):
            trans_df['i->m'].iloc[j-1] = -1000
        
        if (trans_df['d->m'].iloc[j-1] == '*'):
            trans_df['d->m'].iloc[j-1] = -1000
            
        n1 = vm[j-1,i-1] + float(trans_df['m->m'].iloc[j-1])
        n2 = vi[j-1,i-1] + float(trans_df['i->m'].iloc[j-1])
        n3 = vd[j-1,i-1] + float(trans_df['d->m'].iloc[j-1])
        n  = float(match_df[seq[j-1]].iloc[i]) + np.log(seql.count(seq[j-1]))
        maxar = np.array([n1,n2,n3],dtype=np.float64)
        vm[j,i] = n + np.amax(maxar)
        
        # Match trace matrix
        if(np.argmax(maxar) == 0):
            vmt[j,i] = 'd'
        if(np.argmax(maxar) == 1):
            vmt[j,i] = 'u'
        if(np.argmax(maxar) == 2):
            vmt[j,i] = 'l'
        
        # Insert matrix
        if (trans_df['m->i'].iloc[j] == '*'):
            trans_df['m->i'].iloc[j] = -1000
        
        if (trans_df['i->i'].iloc[j] == '*'):
            trans_df['i->i'].iloc[j] = -1000
            
        n1 = vm[j,i-1] + float(trans_df['m->i'].iloc[j])
        n2 = vi[j,i-1] + float(trans_df['i->i'].iloc[j])
#         n3 = vd[j,i-1]
        n  = float(insert_df[seq[j-1]].iloc[i]) + np.log(seql.count(seq[j-1]))
        maxar = np.array([n1,n2],dtype=np.float64)
        vi[j,i] = n + np.amax(maxar)
        
        # Insert trace matrix
        if(np.argmax(maxar) == 0):
            vit[j,i] = 'd'
        if(np.argmax(maxar) == 1):
            vit[j,i] = 'u'
#         if(np.argmax(maxar) == 2):
#             vit[j,i] = 'l'
        
        # Delete matrix
        if (trans_df['m->d'].iloc[j-1] == '*'):
            trans_df['m->d'].iloc[j-1] = -1000
        
        if (trans_df['d->d'].iloc[j-1] == '*'):
            trans_df['d->d'].iloc[j-1] = -1000
            
        n1 = vm[j-1,i] + float(trans_df['m->d'].iloc[j-1])
#         n2 = vi[j-1,i] 
        n3 = vd[j-1,i] + float(trans_df['d->d'].iloc[j-1])
        maxar = np.array([n1,n3],dtype=np.float64)
        vd[j,i] = np.amax(maxar)
        
        # Delete trace matrix
        if(np.argmax(maxar) == 0):
            vdt[j,i] = 'd'
#         if(np.argmax(maxar) == 1):
#             vdt[j,i] = 'u'
        if(np.argmax(maxar) == 1):
            vdt[j,i] = 'l'
        
        # Termination matrix
        maxar = np.array([vm[j,i],vi[j,i],vd[j,i]],dtype=np.float64)
        vt[j,i] = np.amax(maxar)
        
        # Traceback matrix
        if(np.argmax(maxar) == 0):
            tb[j,i] = 'M'
        if(np.argmax(maxar) == 1):
            tb[j,i] = 'I'
        if(np.argmax(maxar) == 2):
            tb[j,i] = 'D'
            
test_df = pd.DataFrame(vdt)
t=test_df.to_csv(test,index=False)

# Trace back
i = tb.shape[0]-1; j = tb.shape[1]-1
st = match_df.columns[1:]

a0 = []; a1 = []; a2 = [];
while(tb[i,j] != 'x'):
    a0.append(seq[i-1])
    maxar = vt[i,:]
    j = np.argmax(maxar)
    a1.append(tb[i,j])
    i = i - 1
    
# Reverse the arrays
as0 = ''.join(a0[::-1]); as1 = ''.join(a1[::-1])
print(as0)
print(as1)
print(as0)