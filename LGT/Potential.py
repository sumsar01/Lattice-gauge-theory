from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
from Momentum import *
import numpy as np 
import os
from Storage import *

# =============================================================================
# Building mass term with gauge field
# =============================================================================
    
def Make_mass(x, y, r, m, a, C):
    
    save_dir = './Data/Operators/Mass/'
    filename = 'Mass_a=' + str(a) + '_m=' + str(m) + '_x=' + str(x) + '_y=' +str(y)

    if os.path.isfile(save_dir + filename + '.p'):
        Mass = load_data(save_dir, filename)
    else:
        Mass = mass_spin(x, y, r, m, a, C)
        save_data(Mass, save_dir, filename)
        
    return Mass

def mass_spin(x, y, r, m, a, C):

    save_dir = './Data/Operators/Simple/'
    filename = 'norm' + '_x=' + str(x) + '_y=' +str(y)
    
    if os.path.isfile(save_dir + filename + '.p'):
        inner_prod = load_data(save_dir, filename)
    else:
        inner_prod = 0
        j = 0
        while j < y:
            i = 0
            while i < x:
                
                inner_prod += ferm_dag_1(x, y, i, j)*ferm_1(x, y, i, j) - ferm_dag_2(x, y, i, j)*ferm_2(x, y, i, j)#norm(x, y, i, j)
                
                i += 1
            j += 1
        save_data(inner_prod, save_dir, filename)
    
    H_m = []
    
    mass = 0
    mass += C[0]*(m)*inner_prod
    mass += C[1]*(m + 2*r/a)*inner_prod
    mass += C[2]*(m + 2*2*r/a)*inner_prod
    H_m.append(mass)
        
    M = sum(H_m)
             
    return M

# =============================================================================
# Building electric potential term
# =============================================================================

def Make_electric_pot(x, y, a, e):
    
    save_dir = './Data/Operators/Electric_potential/'
    filename = 'Electric_term' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_e=' + str(e)
    
    if os.path.isfile(save_dir + filename + '.p'):
        E_pot = load_data(save_dir, filename)
    else:
        E_pot = electric_pot(x, y, a, e)
        save_data(E_pot, save_dir, filename)

    return E_pot


def electric_pot(x, y, a, g):
     
    E_pot = 0
    N = x*y
    J = (-a*g/2)**2       

    save_dir = './Data/Operators/Electric_potential/'
    filename = 'Electric_field' + '_x=' + str(x) + '_y=' + str(y)

    if os.path.isfile(save_dir + filename + '.p'):
        E_pot = load_data(save_dir, filename)
    else:    
        j = 0
        while j < y:
            i = 0
            while i < x:
            
                E_pot += 1/2*Ex(x, y, i, j)*Ex(x, y, i, j)
                E_pot += 1/2*Ey(x, y, i, j)*Ey(x, y, i, j)
            
                i += 1
            j += 1
        save_data(E_pot, save_dir, filename)
    
    E_pot = J*E_pot


    return E_pot

# =============================================================================
# Building magnetic potential term
# =============================================================================

def Make_magnetic_pot(x, y, a, e):
    
    save_dir = './Data/Operators/Magnetic_potential/'
    filename = 'Magnetic_term' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_e=' + str(e)
    
    if os.path.isfile(save_dir + filename + '.p'):
        M_pot = load_data(save_dir, filename)
    else:
        M_pot = magnetic_pot(x, y, a, e)
        save_data(M_pot, save_dir, filename)

    return M_pot

def magnetic_pot(x, y, a, g):
    
    M_pot = 0
    N = x*y
    J = 1/(4*a**4*g**2)

    save_dir = './Data/Operators/Simple/'
    filename = 'magnetic_plaquette' + '_x=' + str(x) + '_y=' + str(y)

    if os.path.isfile(save_dir + filename + '.p'):
        plaquette = load_data(save_dir, filename)
    else:
        j = 0
        while j < y:
            i = 0
            while i < x:
            
                plaquette = Ux(x, y, i, j)*Uy(x, y, i+1, j)*Ux(x, y, i, j+1).dag()*Uy(x, y, i, j).dag()            
                plaquette += (plaquette + plaquette.dag())

                i += 1
            j += 1
        save_data(plaquette, save_dir, filename)


    M_pot = J*plaquette

    return M_pot






















