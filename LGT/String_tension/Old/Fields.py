from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
import numpy as np 
import os
from Storage import *
from Operator_projection import *

#
# Representationen er omvendat af den normale så sigma minus er sigma plus
#
# =============================================================================
# Build gauge fields
# =============================================================================

def U(t, args):
    """ Creates A-field """
    E = args['E']
    a = args['a']
    e = args['e']
    A_field = -E*t
    return np.exp(1.j*e*a*A_field)

def U_dag(t, args):
    """ Creates A-field """
    E = args['E']
    a = args['a']
    e = args['e']
    A_field = -E*t
    return np.exp(-1.j*e*a*A_field)

# Two lattice sites

def U2(t, args):
    """ Creates A-field """
    E = args['E']
    a = args['a']
    e = args['e']
    A_field = -E*t
    return np.exp(1.j*e*a*A_field)*np.exp(1.j*e*a*A_field)

def U2_dag(t, args):
    """ Creates A-field """
    E = args['E']
    a = args['a']
    e = args['e']
    A_field = -E*t
    return np.exp(-1.j*e*a*A_field)*np.exp(-1.j*e*a*A_field)


# =============================================================================
# Build fermionic fields
# =============================================================================

def ferm_dag(x, y, i, j):

    
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N

    
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [sigmam(), sigmam()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)
    
    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_dag_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    
    return field
    
    

def ferm(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [sigmap(), sigmap()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)
    
    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)

    return field
    
# =============================================================================
# Build fermionic fields components
# =============================================================================

def ferm_dag_1(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [sigmam(), qeye(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)
    
    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_dag_1_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    return field
    
def ferm_dag_2(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
   
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [qeye(2), sigmam()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)
    
    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_dag_2_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)
    
    
    return field
    

def ferm_1(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [sigmap(), qeye(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)
    
    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_1_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)

    return field

def ferm_2(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N    
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [qeye(2), sigmap()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = field + void(x, y)

    save_dir = './Data/Operators/Simple/'
    filename = 'ferm_2_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
    
    if os.path.isfile(save_dir + filename + '.p'):
        field = load_data(save_dir, filename)
    else:
        field = project_op(x, y, field)
        save_data(field, save_dir, filename)

    return field


"""

def norm(x, y, i, j):

    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    
    if N < 1:
        raise ValueError("integer N must be larger or equal to 1")

    U = [- tensor(sigmaz(), identity(2))*1/2 + tensor(identity(2), sigmaz())*1/2]

    return tensor([identity(2)] * (target) * 2 + U +
                  [identity(2)] * (N - target - 1) * 2)

"""





















