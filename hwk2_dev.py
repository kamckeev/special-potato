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

lines = [line.strip() for line in lines]
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
            #fwd_matrix[0,1] = 0.0000000000000000000000001
            fwd_matrix[0,1] = -1000000000000000000000000
        else:
            prev_AT = fwd_matrix[i-1,0] #Previous probability in AT region
            prev_GC = fwd_matrix[i-1,1] #Previous probability in GC region
            max_p = max(prev_AT, prev_GC)

            #A/T in AT trans_p prob 
            trans_p00 = math.log((emission_p[0,base])*(trans_p[0,0])) 
            
            #G/C in AT trans_p prob
            trans_p10 = math.log((emission_p[0,base])*(trans_p[1,0])) 
            
            #A/T in GC trans_p prob
            trans_p01 = math.log((emission_p[1,base])*(trans_p[0,1])) 
            
            #G/C in GC trans_p prob
            trans_p11 = math.log((emission_p[1,base])*(trans_p[1,1])) 
            
            fwd_matrix[i,0]=max_p+math.log(math.exp(-max_p+prev_AT
                                                      +trans_p00)
                                             +math.exp(-max_p 
                                                       + prev_GC+trans_p10))
            fwd_matrix[i,1]=max_p+math.log(math.exp(-max_p
                                                      +prev_AT+trans_p01)
                                             +math.exp(-max_p 
                                                       + prev_GC+trans_p11))
            
            final_AT = fwd_matrix[lenseq-1, 0] 
            final_GC= fwd_matrix[lenseq-1, 1]
    final_p=final_AT+math.log(1+(math.exp(final_GC-final_AT)))
    #return(f"The log-likelihood is {final_p}")
    return(final_p)
#%%
"""
A sequence of length 1200 was generated with a hidden Markov model that allowed
some regions to be GC-rich and some regions to be  AT-rich.  The model used to 
generate the data was very similar to the one that we studied in class...

Position 1 in the sequence was forced to be in a AT-rich region.

Thereafter, the probability that position i+1 was in a GC-rich region given 
that position i was in a GC-rich region was set to 0.95. The probability that
position i+1 was in an AT-rich region given that position i was in an AT-rich
region was set to 0.98.

In GC-rich regions, the probabilities of G and C occupying sites were each 0.3 
whereas the probabilities of A and T in these GC-rich regions were each 0.2.

In AT-rich regions, the probabilities of A and T occupying sites were each 0.3 
whereas the probabilities of G and C in these AT-rich regions were each 0.2.
"""
beg_state=np.array((1,0))
emission_p=np.array(((0.3,0.3,0.2,0.2), (0.2,0.2,0.3,0.3)))
trans_p=np.array(((0.98, 0.02), (0.05, 0.95)))
#%%
"""
1.  Implement the forward algorithm to calculate the logarithm of the 
probability of observing this sequence given the model.  When calculating this 
log-likelihood, do it for the true values of the parameters that are listed above.
What log-likelihood value results?
"""

#print(trans_p)
print(f"#1: {fwd_alg(sequence, trans_p, emission_p, beg_state)}")
#%%
"""
2.  Implement the forward algorithm again to analyze the same data set.  
But, this time use the true values of the parameters except 0.5 should be the
probability that position i+1 was in a GC-rich region given that position i 
was in a GC-rich region and 0.8 should be the probability that position i+1 was
in an AT-rich region given that position i was in an AT-rich region.  What 
value of the log-likelihood results?
"""
trans_p=np.array(((0.8, 0.2), (0.5, 0.5)))
print(f"#2: {fwd_alg(sequence, trans_p, emission_p, beg_state)}")
#%%
"""
3. Implement the forward algorithm again to analyze the same data set.  
But, this time use the true values of the parameters except 0.51 should be the 
probability that position i+1 was in a GC-rich region given that position i 
was in a GC-rich region and 0.51 should be the probability that position i+1 
was in an AT-rich region given that position i was in an AT-rich region.  
What value of the log-likelihood results?
"""
trans_p=np.array(((0.51, 0.49), (0.51, 0.49)))
print(f"#3: {fwd_alg(sequence, trans_p, emission_p, beg_state)}")
#%%
def viterbi_alg(sequence, trans_p, emision_p, beg_state):
    v_matrix = np.zeros((lenseq, states)) #initialize prob matrix with zeros
    path_matrix = np.zeros((lenseq, states)) #initialize path matrix with zeros
    for i in range(lenseq):
        if sequence[i] == 'A':
            base = 0
        elif sequence[i] == 'T':
            base = 1
        elif sequence[i] == 'G':
            base = 2
        elif sequence[i] == 'C':
            base = 3
            
        v_matrix[0,0] = math.log((beg_state[0]*emission_p[0,base])) #prob =(inital)*(emission)
        #trying negative number
        v_matrix[0,1] = -1000000000000000000000000000000000
        #trying really small number
        #v_matrix[0,1] = 0.0000000000000000000000000000000001
        
        prev_AT = v_matrix[i-1,0]
        prev_GC = v_matrix[i-1,1]

        trans00 = prev_AT+math.log(emission_p[0,base])+math.log(trans_p[0,0]) 
        trans10 = prev_GC+math.log(emission_p[0,base])+math.log(trans_p[1,0]) 
        trans01 = prev_AT+math.log(emission_p[1,base])+math.log(trans_p[0,1])
        trans11 = prev_GC+math.log(emission_p[1,base])+math.log(trans_p[1,1]) 

        #Switch to GC
        if trans00 < trans10: 
            path_matrix[i,0] = 1 
            v_matrix[i,0] = trans10 
            
        #Stay in AT
        else: 
            path_matrix[i,0] = 0 
            v_matrix[i,0] = trans00 
            
        #Switch to AT    
        if trans01 > trans11: 
            path_matrix[i,1] = 0 
            v_matrix[i,1] = trans01 
            
        else: #Stay in GC
            path_matrix[i,1] = 1 
            v_matrix[i,1] = trans11  
            
    if v_matrix[lenseq-1,0] > v_matrix[lenseq-1,1]:
        max_p = (v_matrix[lenseq-1,0]) 
    else:
        max_p = (v_matrix[lenseq-1,1]) 
    
    #Ending AT rich
    if v_matrix[lenseq-1,0] > v_matrix[lenseq-1,1]:
        state = 0  
    #Ending GC rich
    else:
        state = 1 
    

#traceback pathway
    final_path = [state] #Start with whatever the last state was
    for j in range (lenseq-2): 
        add = int(path_matrix[lenseq-j-2, state]) 
        final_path.insert(0, add) 
        state = int(path_matrix[lenseq-j-2, state])

    return f"The resulting log-likelihood value is {max_p}\r \
This is the final path:\n {str(final_path)}" 
#%%
"""
4.  For the parameters values listed in Question 1, implement the Viterbi 
algorithm to analyze the simulated sequence data.  The implementation should 
report the most probable path of AT-rich and GC-rich states through the hidden
Markov model.  It should also report the logarithm of the probability of 
generating the sequence data with this specific path.
"""
trans_p=np.array(((0.98, 0.02), (0.05, 0.95)))
print("For the following questions, 0 indicates an AT rich state and \
      1 indicates a GC rich state\n")
print(f"#5: {viterbi_alg(sequence, trans_p, emission_p, beg_state)}")
#%%
"""
5.  Answer Question 4 again except now use the parameter values that were used 
in Question 2.  In a few sentences, try to explain why the reconstructed path 
from this question differs from the path reconstructed for Question 4.
"""
trans_p=np.array(((0.8, 0.2), (0.5, 0.5)))
print(f"#5: {viterbi_alg(sequence, trans_p, emission_p, beg_state)}\r")
print("This reconstructed path is different from the path in #4 \
because we used different tranistion probabilities.\nnew sentence.\
          \nnew sentence")
#%%
"""
6.  Answer Question 4 again except now use the parameter values that were used 
in Question 3.  In a few sentences, try to explain why the reconstructed path 
from this question differs from the path reconstructed for Question 4.
"""
trans_p=np.array(((0.51, 0.49), (0.51, 0.49)))
print(f"#6: {viterbi_alg(sequence, trans_p, emission_p, beg_state)}\r")
print("This reconstructed path is different from the path in #4 \
because we used different tranistion probabilities.\nnew sentence.\
          \nnew sentence")





