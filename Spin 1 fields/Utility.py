from qutip import *
import numpy as np 
import os
from Storage import *

###############################################################################
# Operator product
###############################################################################

def ope_prod(ope_1, ope_2, ope_3, ope_4):
    
    length = len(ope_1)
    i = 0
    ope = []
    
    if ope_4 != 0:
        while i < length:
            ope.append(ope_1[i]*ope_2[i]*ope_3[i]*ope_4[i])
            i += 1
            
    if ope_3 != 0 and ope_4 == 0:
        while i < length:
            ope.append(ope_1[i]*ope_2[i]*ope_3[i])
            i += 1

    if ope_3 == 0 and ope_4 == 0:
        while i < length:
            ope.append(ope_1[i]*ope_2[i])
            i += 1        
            
    return ope

def find_ternary(num):  #2
    if num == 0:
        return "0"
    quotient = num/3    #3
    remainder = num%3
    if quotient == 0:   #4
        return ""
    else:
        return find_ternary(int(quotient)) + str(int(remainder))    #5
    number = int(input("Enter a number : ")) #1
    return int(find_ternary(number))






























