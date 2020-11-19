from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
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
    filename = 'Wilson_loop_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)
    
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

# disorder operator

def Hooft_string(x, y, i, j, phase):

    save_dir = './Data/Operators/Hooft/'
    filename = 'Hooft_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)

    if os.path.isfile(save_dir + filename + '.p'):
        Y = load_data(save_dir, filename)
    else:
        E1 = project_op(x, y, Ex(x, y, i, j))
        E2 = project_op(x, y, Ex(x, y, i+1, j))
        Y = (1.j*phase*E1).expm()*(1.j*phase*E2).expm()
        save_data(Y, save_dir, filename)
    
    return Y

def Wilson_line(x, y, i, j):

    save_dir = './Data/Operators/Wilson_line/'
    filename = 'line' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        line = load_data(save_dir, filename)
    else:
        line = ope_prod(ferm_dag_1(x, y, i, j), Ux(x, y, i, j), ferm_2(x, y, i+1, j), 0)
        line = project_op(x, y, line)
        save_data(line, save_dir, filename)

    return line
















