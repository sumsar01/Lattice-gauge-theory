from qutip import *
from itertools import product
import itertools as itertools
from pathlib import Path
import numpy as np 
import time
from Utility import *

def advanced_projector(x, y):
    
    Nm = (x*y)*2
    G_n = [0]*x*y
    good_stateReps = []
    
    
    def Generator_on_site(Rep):
        charge = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                P = (j*x + i)*2
                A = P+1
                i1 = (i-1) % x
                j1 = (j-1) % y
                E_x = P
                E_y = A
                E_x_back = (j*x + i1)*2
                E_y_back = (j1*x + i)*2 
                charge.extend([(Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
                i += 1
            j += 1
        return charge
    
    starttime = time.time()
    j = 0
    while j < 3**(Nm):
        num = int(find_ternary(j))
        stateRep_gauge = [int(x)-1 for x in '{0:0{1}}'.format(num,Nm)]
            
        stateRep = stateRep_gauge
              
        if Generator_on_site(stateRep) == [0]*x*y:
            good_stateReps.append(stateRep)
        
        progress = round(j/2**(Nm),4)
        speed = round(progress/(time.time() - starttime),6)
        print('\rN = ' + str(Nm) + ', progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='') 
            
        j += 1
  
    
    return good_stateReps 
      
   












          
                
                
                
                