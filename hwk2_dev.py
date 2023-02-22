# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 15:37:44 2023

@author: jmcco
"""
import math
import numpy as np
#%%
""" Importing sequence files"""
import_file=open('hwk2_sequence.txt', 'r')
lines=import_file.readlines()

#lines = [line.strip() for line in lines]

for line in lines:
    line = line.strip()
sequence = ''.join(lines)

import_file.close()
#%%
#Forward Algorithm
lenseq=len(sequence)
states=2
def fwd_alg(sequence, trans_p, emission_p, beg_state): 
    #Initialize probability matrix with zeros
    fwd_matrix = np.zeros((lenseq, states)) 
    for i in range(lenseq):   
        if sequence[i] == 'A': 
            base = 0
        elif sequence[i] == 'T':
            base = 1
        elif sequence[i] == 'G':
            base = 2
        elif sequence[i] == 'C':
            base = 3    
        if i == 0: #Use iniital probabilities
            fwd_matrix[0,0] = math.log((beg_state[0]*emission_p[0,base]))
            fwd_matrix[0,1] = -1000000000000000000
        else:
            prev_AT = fwd_matrix[i-1,0] #Previous probability in AT region
            prev_GC = fwd_matrix[i-1,1] #Previous probability in GC region
            maxProb = max(prev_AT, prev_GC)

            #A/T in AT trans_p prob 
            trans_p00 = math.log((emission_p[0,base])*(trans_p[0,0])) 
            
            #G/C in AT trans_p prob
            trans_p10 = math.log((emission_p[0,base])*(trans_p[1,0])) 
            
            #A/T in GC trans_p prob
            trans_p01 = math.log((emission_p[1,base])*(trans_p[0,1])) 
            
            #G/C in GC trans_p prob
            trans_p11 = math.log((emission_p[1,base])*(trans_p[1,1])) 
            
            fwd_matrix[i,0]=maxProb+math.log(math.exp(-maxProb+prev_AT
                                                      +trans_p00)
                                             +math.exp(-maxProb 
                                                       + prev_GC+trans_p10))
            fwd_matrix[i,1]=maxProb+math.log(math.exp(-maxProb
                                                      +prev_AT+trans_p01)
                                             +math.exp(-maxProb 
                                                       + prev_GC+trans_p11))
            
            final_AT = fwd_matrix[lenseq-1, 0] 
            final_GC= fwd_matrix[lenseq-1, 1]
    final_prob=final_AT+math.log(1+(math.exp(final_GC-final_AT)))
    print("The log-likelihood is",final_prob)
#%%
#1
beg_state=np.array((1,0))
emission_p=np.array(((0.3,0.3,0.2,0.2), (0.2,0.2,0.3,0.3)))
trans_p=np.array(((0.98, 0.02), (0.05, 0.95)))
print(trans_p)
print("#1")
fwd_alg(sequence, trans_p, emission_p, beg_state)
#%%
#2
trans_p=np.array(((0.8, 0.2), (0.5, 0.5)))
print("#2")
fwd_alg(sequence, trans_p, emission_p, beg_state)
#%%
#3
trans_p=np.array(((0.51, 0.49), (0.51, 0.49)))
print("#3")
fwd_alg(sequence, trans_p, emission_p, beg_state)
