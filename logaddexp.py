# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 17:53:35 2023

@author: jmcco
"""
import numpy as np
import math


prob1 = np.log(1e-50)
prob2 = np.log(2.5e-50)
prob12 = np.logaddexp(prob1, prob2)
prob12

print(np.exp(prob12))
print(np.exp(prob1 + prob2))
#%%
max_p=A=1
prev_AT=B=1
trans_p00=C=0.2
prev_GC=D= 5
trans_p10=E=0.8
trans_p01=F=0.5
trans_p11=G=0.5
final_AT=H=2
final_GC=I=1
#%%
#fwd_matrix[i,0]=max_p+math.log(math.exp(-max_p+prev_AT+trans_p00) + math.exp(-max_p + prev_GC+trans_p10))
orig=max_p+math.log(math.exp(-max_p+prev_AT+trans_p00) + math.exp(-max_p + prev_GC+trans_p10))
orig2=max_p+math.log(math.exp(-A+B+C) + math.exp(-A + D +E))
print(orig,orig2)
#%%
orig2=max_p+math.log(math.exp(-A+B+C) + math.exp(-A + D +E))
#orig4=np.logaddexp(max_p,(math.exp(-A+B+C) + math.exp(-A + D +E)))
#print(orig2,orig4)
#orig3=max_p+math.log(math.exp(np.logaddexp.reduce([-A,B,C]))+math.exp(np.logaddexp.reduce([-A,D,E])))
#print(orig2,orig3)
#%%
pt1=max_p+math.log(math.exp(-max_p+prev_AT+trans_p00) + math.exp(-max_p + prev_GC+trans_p10))
interim1 = math.exp(B+C-A)
interim2 = math.exp(D+E-A)
eff = max_p+math.log(np.logaddexp(interim1, interim2))

eff2 = max_p+np.logaddexp((B+C-A),(D+E-A))

print(interim1)
print(interim2)
print(eff2)
print(pt1)
#%%
pt2=max_p+math.log(math.exp(-max_p+prev_AT+trans_p01)+math.exp(-max_p + prev_GC+trans_p11))
interim3=math.exp(+B+F-A)
interim4=math.exp(D+G-A)
gee=max_p+math.log(np.logaddexp(interim3,interim4))

print(interim3)
print(interim4)
print(gee)
print(pt2)
#%%
pt3=final_AT+math.log(1+(math.exp(final_GC-final_AT)))
interim5=math.exp(I-H)
ayy=H+math.log(1+interim5)
print(pt3)
print(interim5)
print(ayy)
#%%
final1 = np.logaddexp(max_p,np.logaddexp(interim1,interim2))
print(final1)
print(orig)
#%%
prev_AT=J=1
prev_GC=K=2
#%%
L=math.log(emission_p[0,base])
M=math.log(emission_p[1,base])
N=math.log(trans_p[0,0])
O=math.log(trans_p[1,0])
P=math.log(trans_p[0,1])
Q=math.log(trans_p[1,1])
#%%
trans00 = prev_AT+math.log(emission_p[0,base])+math.log(trans_p[0,0]) 
trans10 = prev_GC+math.log(emission_p[0,base])+math.log(trans_p[1,0]) 
trans01 = prev_AT+math.log(emission_p[1,base])+math.log(trans_p[0,1])
trans11 = prev_GC+math.log(emission_p[1,base])+math.log(trans_p[1,1])
#%%
#test1 = prev_AT+math.log(emission_p[0,base])+math.log(trans_p[0,0]) 
test2=prev_AT+
