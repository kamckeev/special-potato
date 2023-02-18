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