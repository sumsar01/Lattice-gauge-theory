from Phase import *
from qutip import *
from Fields import *
from Fields_spin1 import *
import numpy as np 
import os
from Storage import *
from Utility import *

# =============================================================================
# Building mass term with gauge field
# =============================================================================

def mass_spin(x, y, r, m, a, C):

    save_dir = './Data/Operators/Mass/non_periodic/'
    filename = 'Mass_a=' + str(a) + '_m=' + str(m) + '_x=' + str(x) + '_y=' +str(y)
 
    if os.path.isfile(save_dir + filename + '.p'):
        H_mass = load_data(save_dir, filename)
    else:
        Mass = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                save_dir2 = './Data/Operators/Simple/non_periodic/'
                filename_P = 'Site_mass_P=' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
                filename_A = 'Site_mass_A=' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
                target = (j*y + i)*2
                
                if os.path.isfile(save_dir2 + filename_P + '.p'):
                    P = load_data(save_dir2, filename_P)
                else:
                    P = ope_prod(ferm_dag_1(x, y, i, j), ferm_1(x, y, i, j), 0, 0)
                    P = project_op(x, y, P)
                    save_data(P, save_dir2, filename_P)
                    
                if os.path.isfile(save_dir2 + filename_A + '.p'):
                    A = load_data(save_dir2, filename_A)
                else:
                    A = ope_prod(ferm_dag_2(x, y, i, j), ferm_2(x, y, i, j), 0, 0)
                    A[target+1] = -A[target+1]
                    A = project_op(x, y, A)
                    save_data(A, save_dir2, filename_A)
                
                term = P + A
                Mass.append(term)
                
                i += 1
            j += 1
        H_mass = sum(Mass)
        H_mass = H_mass*C[1]*(m + 2*r/a)
        save_data(H_mass, save_dir, filename)
             
    return H_mass

# =============================================================================
# Building electric potential term
# =============================================================================

def electric_pot(x, y, a, g):

    save_dir = './Data/Operators/Electric_potential/non_periodic/'
    filename = 'Electric_term' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_e=' + str(g)     
    E_pot = 0
    N = x*y
    J = 1/2*(-a*g/2)**2       

    if os.path.isfile(save_dir + filename + '.p'):
        H_E = load_data(save_dir, filename)
    else:
        save_dir2 = './Data/Operators/Simple/non_periodic/'
        filename2 = 'lattice_Sz_' + '_x=' + str(x) + '_y=' + str(y)

        if os.path.isfile(save_dir2 + filename2 + '.p'):
                    H_E = load_data(save_dir2, filename2)
        else:
            E_field = []
            j = 0
            i = 0
            E_pot_x = ope_prod(Ex(x, y, i, j), Ex(x, y, i, j), 0, 0)
            E_pot_x = project_op(x, y, E_pot_x)

            E_pot_y = ope_prod(Ey(x, y, i, j), Ey(x, y, i, j), 0, 0)
            E_pot_y = project_op(x, y, E_pot_y)
            term = E_pot_x + E_pot_y
            E_field.append(term)
            
            j = 0
            i = 1
            E_pot_y = ope_prod(Ey(x, y, i, j), Ey(x, y, i, j), 0, 0)
            E_pot_y = project_op(x, y, E_pot_y)
            term = E_pot_y
            E_field.append(term)
        
            j = 1
            i = 0
            E_pot_x = ope_prod(Ex(x, y, i, j), Ex(x, y, i, j), 0, 0)
            E_pot_x = project_op(x, y, E_pot_x)
            term = E_pot_x
            E_field.append(term)
            H_E = sum(E_field)
            save_data(H_E, save_dir2, filename2)
            
        H_E = J*H_E
        save_data(H_E, save_dir, filename)

    return H_E

# =============================================================================
# Building magnetic potential term
# =============================================================================



def magnetic_pot(x, y, a, g):
 
    save_dir = './Data/Operators/Magnetic_potential/non_periodic/'
    filename = 'Magnetic_term' + '_x=' + str(x) + '_y=' + str(y) + '_a=' + str(a) + '_e=' + str(g)
    M_pot = 0
    N = x*y
    J = 1/(4*a**4*g**2)

    if os.path.isfile(save_dir + filename + '.p'):
        H_B = load_data(save_dir, filename)
    else:
        B_field = []
        j = 0
        i = 0
        save_dir2 = './Data/Operators/Simple/non_periodic/'
        filename_R = 'plaquette_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)
        filename_L = 'plaquette_dag_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)

        if os.path.isfile(save_dir2 + filename_R + '.p'):
            plaquette = load_data(save_dir2, filename_R)
        else:
            plaquette = ope_prod(Ux(x, y, i, j), Uy(x, y, i+1, j), Ux_dag(x, y, i, j+1), Uy_dag(x, y, i, j))
            plaquette = project_op(x, y, plaquette)
            save_data(plaquette, save_dir2, filename_R)                    

        term = (plaquette + plaquette.dag())
        B_field.append(term)

        H_B = sum(B_field)
        H_B = J*H_B
        save_data(H_B, save_dir, filename)

    return H_B

