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
                P = string_op(x, y, P)
                A = string_op(x, y, A)
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
        P = string_op(x, y, P)
        A = string_op(x, y, A)        
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
        loop = string_op(x, y, loop)
        save_data(loop, save_dir, filename)

    return loop

def Wilson_loop2(x, y, i, j):
  
    save_dir = './Data/Operators/Wilson_loop/'
    filename = 'Wilson_loop2_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        loop = load_data(save_dir, filename)
    else:
        loop = ope_prod(Ux(x, y, i, j), Ux(x, y, i+1, j), 0, 0)
        loop = string_op(x, y, loop)
        save_data(loop, save_dir, filename)

    return loop

def Wilson_loop3(x, y, i, j):

    save_dir = './Data/Operators/Wilson_loop/'
    filename = 'Wilson_loop2_x=' + str(x) + '_y=' +str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        loop = load_data(save_dir, filename)
    else:
        loop = ope_prod(Uy(x, y, i, j), Uy(x, y, i, j+1), 0, 0)
        loop = string_op(x, y, loop)
        save_data(loop, save_dir, filename)
    
    return loop




def electric(x, y, a):

    save_dir = './Data/Operators/Electric_potential/'
    filename = 'Electric' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a)     
    E_pot = 0
    N = x*y

    if os.path.isfile(save_dir + filename + '.p'):
        H_E = load_data(save_dir, filename)
    else:
        E_field = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                save_dir2 = './Data/Operators/Simple/'
                filename_x = 'Link_Sz_x_single_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
                filename_y = 'Link_Sz_y_single_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
                
                if os.path.isfile(save_dir2 + filename_x + '.p'):
                    E_pot_x = load_data(save_dir2, filename_x)
                else:
                    E_pot_x = string_op(x, y, Ex(x, y, i, j))
                    save_data(E_pot_x, save_dir2, filename_x)
               
                if os.path.isfile(save_dir2 + filename_y + '.p'):
                    E_pot_y = load_data(save_dir2, filename_y)
                else:
                    E_pot_y = string_op(x, y, Ey(x, y, i, j))
                    save_data(E_pot_y, save_dir2, filename_y)
                
                term = E_pot_x + E_pot_y
                E_field.append(term)
            
                i += 1
            j += 1
        H_E = sum(E_field)
        save_data(H_E, save_dir, filename)

    return H_E


def Wilson_line(x, y, i, j):

    save_dir = './Data/Operators/Wilson_line/'
    filename = 'line' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        line = load_data(save_dir, filename)
    else:
        line = ope_prod(ferm_dag_1(x, y, 0, 0), ferm_2(x, y, 2, 1), Ux(x, y, 0, 0), Ux(x, y, 1, 0))
        line = ope_prod(line, Uy(x, y, 2, 0), 0, 0)
        line = string_op(x, y, line)
        save_data(line, save_dir, filename)

    return line



def electric_x(x, y, a, i, j):

    save_dir = './Data/Operators/Electric_potential/'
    filename = 'Electric_x' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_i=' + str(i) + '_j=' + str(j)     
    E_pot = 0
    N = x*y

    if os.path.isfile(save_dir + filename + '.p'):
        H_E = load_data(save_dir, filename)
    else:
        H_E = string_op(x, y, Ex(x, y, i, j))
        save_data(H_E, save_dir, filename)

    return H_E

def electric_y(x, y, a, i, j):

    save_dir = './Data/Operators/Electric_potential/'
    filename = 'Electric_y' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_i=' + str(i) + '_j=' + str(j)     
    E_pot = 0
    N = x*y

    if os.path.isfile(save_dir + filename + '.p'):
        H_E = load_data(save_dir, filename)
    else:
        H_E = string_op(x, y, Ey(x, y, i, j))
        save_data(H_E, save_dir, filename)

    return H_E












