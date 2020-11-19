from Phase import *
from qutip import *
from Fields import *
from Fields_spin1 import *
from Momentum_spin import *
from Potential import *
import numpy as np 

# =============================================================================
# Number operator
# =============================================================================

def num_ope(x, y):
    
    save_dir = './Data/Operators/Number/'
    filename = 'number_ope_x=' + str(x) + '_y=' +str(y)
    
    if os.path.isfile(save_dir + filename + '.p'):
        number_ope = load_data(save_dir, filename)
    else:
        number_ope = [] 
        j = 0
        while j < y:
            i = 0
            while i < x:
                P = ope_prod(ferm_dag_1(x, y, i, j), ferm_1(x, y, i, j), 0, 0)
                A = ope_prod(ferm_2(x, y, i, j), ferm_dag_2(x, y, i, j), 0, 0)
                P = project_op(x, y, P)
                A = project_op(x, y, A)
                num = P + A
                number_ope.append(num)
           
                i += 1
            j += 1
        number_ope = sum(number_ope)
        save_data(number_ope, save_dir, filename)

    return number_ope

def num_ope_site(x, y, i, j):
    
    save_dir = './Data/Operators/Number/'
    filename = 'number_ope_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)

    if os.path.isfile(save_dir + filename + '.p'):
        number_ope = load_data(save_dir, filename)
    else:
        P = ope_prod(ferm_dag_1(x, y, i, j), ferm_1(x, y, i, j), 0, 0)
        A = ope_prod(ferm_2(x, y, i, j), ferm_dag_2(x, y, i, j), 0, 0)
        P = project_op(x, y, P)
        A = project_op(x, y, A)        
        number_ope = P + A
        save_data(number_ope, save_dir, filename)

    return number_ope



# Loops

def Wilson_loop(x, y, i, j):

    save_dir = './Data/Operators/Wilson_loop/'
    filename = 'plaquette_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        loop = load_data(save_dir, filename)
    else:
        loop = ope_prod(Ux(x, y, i, j), Uy(x, y, i+1, j), Ux_dag(x, y, i, j+1), Uy_dag(x, y, i, j))
        loop = project_op(x, y, loop)
        save_data(loop, save_dir, filename)

    return loop

def Wilson_loop2(x, y, i, j):
  
    save_dir = './Data/Operators/Wilson_loop/'
    filename = 'Wilson_loop2_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        loop = load_data(save_dir, filename)
    else:
        loop = ope_prod(Ux(x, y, i, j), Ux(x, y, i+1, j), 0, 0)
        loop = project_op(x, y, loop)
        save_data(loop, save_dir, filename)

    return loop

def Wilson_loop3(x, y, i, j):

    save_dir = './Data/Operators/Wilson_loop/'
    filename = 'Wilson_loop2_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        loop = load_data(save_dir, filename)
    else:
        loop = ope_prod(Uy(x, y, i, j), Uy(x, y, i, j+1), 0, 0)
        loop = project_op(x, y, loop)
        save_data(loop, save_dir, filename)
    
    return loop


























