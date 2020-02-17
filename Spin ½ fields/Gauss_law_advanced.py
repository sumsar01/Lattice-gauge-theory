from qutip import *
from itertools import product
import itertools as itertools
from pathlib import Path
import numpy as np 
import time

def advanced_projector(x, y):
    
    Nm = (x*y)*4
    G_n = [0]*x*y
    good_stateReps = []
    
    
    def Generator_on_site(Rep):
        charge = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                P = (j*y + i)*2
                A = P+1
                N_tot = x*y*2
                i1 = (i-1) % x
                j1 = (j-1) % y
                E_x = P + N_tot
                E_y = A + N_tot
                E_x_back = (j*y + i1)*2 + N_tot
                E_y_back = (j1*y + i)*2 + N_tot +1
                charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
                i += 1
            j += 1
        return charge
    
    i = 0
    starttime = time.time()
    while i < 2**(Nm):
        stateRep = [int(x) for x in '{0:0{1}b}'.format(i,Nm)]
        
        i += 1
        if i%1000 == 0:
            progress = round(i/2**(Nm)*100,4)
            speed = round(progress/(time.time() - starttime),6)
            print('\rN = ' + str(Nm) + ', progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='')        
        
        
        if Generator_on_site(stateRep) == [0]*x*y:
            good_stateReps.append(stateRep)
    
    return good_stateReps 
      
   
def advanced_projector_number(x, y, n):
    
    Nm = (x*y)*4
    G_n = [0]*x*y
    good_stateReps = []
    
    
    def Generator_on_site(Rep, n):
        charge = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                P = (j*y + i)*2
                A = P+1
                N_tot = x*y*2
                i1 = (i-1) % x
                j1 = (j-1) % y
                E_x = P + N_tot
                E_y = A + N_tot
                E_x_back = (j*y + i1)*2 + N_tot
                E_y_back = (j1*y + i)*2 + N_tot +1
                charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
                if Rep[N] == 1:
                    number += 1
                if Rep[N+1] == 0:
                    number += 1
                    
                i += 1
            j += 1
        if number != n:
            charge = [1, 1, 1, 1]
        return charge
    
    i = 0
    starttime = time.time()
    while i < 2**(Nm):
        stateRep = [int(x) for x in '{0:0{1}b}'.format(i,Nm)]
        
        i += 1
        if i%1000 == 0:
            progress = round(i/2**(Nm)*100,4)
            speed = round(progress/(time.time() - starttime),6)
            print('\rN = ' + str(Nm) + ', progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='')        
        
        
        if Generator_on_site(stateRep) == [0]*x*y:
            good_stateReps.append(stateRep)
                
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Lattice_Projection/'
    filename = 'good_stateReps_N=' + str(Nm) + '_x=' + str(x) + '_y=' +str(y)
    np.save(save_dir + filename,good_stateReps)
    
    return good_stateReps 












          
                
                
                
                