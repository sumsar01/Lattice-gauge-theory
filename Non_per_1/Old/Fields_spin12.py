from Phase import *
from qutip import *
from Fields import *
import numpy as np 
import os
from Storage import *
from Operator_projection import *

# =============================================================================
# Build gauge fields (Spin-1)
# =============================================================================

def Ux(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [sigmam(), identity(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field
    
    save_dir = './Data/Operators/Simple/'
    filename = 'Ux_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    
    return field

def Uy(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(2), sigmam()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field
    
    save_dir = './Data/Operators/Simple/'
    filename = 'Uy_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)    
    
    return field

def void(x, y):
    
    N = x*y
    
    return ([identity(2)] * N * 2)
    
def Ex(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    U = [sigmaz(), identity(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field
    
    save_dir = './Data/Operators/Simple/'
    filename = 'Ex_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    return field
    
def Ey(x, y, i, j):    
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
        
    N = x*y
    target = (j*y + i) % N
    
    U = [identity(2), sigmaz()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field
    
    save_dir = './Data/Operators/Simple/'
    filename = 'Ey_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    return field
    
    
    
    
    
    
    
    
    
    
    
