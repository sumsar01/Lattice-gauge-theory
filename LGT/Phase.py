from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
import numpy as np 

# =============================================================================
# Build hopping phase
# =============================================================================
 

# =============================================================================
# Nearest neighboor phase
# =============================================================================

def x_phase(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    row = [identity(2)] * x * y * 2

    i1 = i + 1 % x
    
    f = j - 2 % y
    while f < y:
        row[(f*x + i) * 2] = row[(f*x + i) * 2]*(-sigmaz())
        row[(f*x + i) * 2 + 1] = row[(f*x + i) * 2 + 1]*(-sigmaz())
        f += 1    
    
    f = 0
    while f < j-1 % y:
        row[(f*x + i1) * 2] = row[(f*x + i1) * 2]*(-sigmaz())
        row[(f*x + i1) * 2 + 1] = row[(f*x + i1) * 2 + 1]*(-sigmaz())
        f += 1
        

        
    phase = row + void(x, y)
    
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Operators/Simple/'
    filename = 'x_phase_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        phase = load_data(save_dir, filename)
    else:
        phase = project_op(x, y, phase)
        save_data(phase, save_dir, filename)  
        
    return phase



def y_phase(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = j*y + i
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [-sigmaz(), -sigmaz()]

    phase = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2) + void(x, y)
    
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Operators/Simple/'
    filename = 'y_phase_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        phase = load_data(save_dir, filename)
    else:
        phase = project_op(x, y, phase)
        save_data(phase, save_dir, filename)  
    
    return phase


# =============================================================================
# Next-nearest neighboor phase
# =============================================================================

def x_phase2(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    row = [identity(2)] * x * y * 2

    i1 = i + 1 % x
    i2 = i + 2 % x
    
    f = j - 2 % y
    while f < y:
        row[(f*x + i) * 2] = row[(f*x + i) * 2]*(-sigmaz())
        row[(f*x + i) * 2 + 1] = row[(f*x + i) * 2 + 1]*(-sigmaz())
        f += 1

    f = 0
    f1 = (f*x + i1) % x*y
    while f < y:
        row[(f1) * 2] = row[(f1) * 2]*(-sigmaz())
        row[(f1) * 2 + 1] = row[(f1) * 2 + 1]*(-sigmaz())
        f += 1

    
    f = 0
    f2 = (f*x + i2) % x*y
    while f < j-1 % y:
        row[(f*x + i2) * 2] = row[(f*x + i2) * 2]*(-sigmaz())
        row[(f*x + i2) * 2 + 1] = row[(f*x + i2) * 2 + 1]*(-sigmaz())
        f += 1
        

        
    phase = row + void(x, y)
    
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Operators/Simple/'
    filename = 'x_phase2_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        phase = load_data(save_dir, filename)
    else:
        phase = project_op(x, y, phase)
        save_data(phase, save_dir, filename)  
        
    return phase


def y_phase2(x, y, i, j):
    
    j2 = j + 1
    
    if j >= y:
        j = j % y
    
    if j2 >= y:
        j2 = j2 % y
    
    N = x*y
    target1 = (j*y+i) * 2
    target2 = (j2*y + i) * 2
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = -sigmaz()
    
    row = [identity(2)] * N *2
    row[target1] = U
    row[target1 + 1] = U
    row[target2] = U
    row[target2 + 1] = U
    
    phase = row + void(x, y)
    
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Operators/Simple/'
    filename = 'y_phase2_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        phase = load_data(save_dir, filename)
    else:
        phase = project_op(x, y, phase)
        save_data(phase, save_dir, filename)  
    
    return phase


















