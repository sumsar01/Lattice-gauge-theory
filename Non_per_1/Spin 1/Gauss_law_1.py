from qutip import *
from itertools import product
import itertools as itertools
import numpy as np 




# =============================================================================
# Building Gauss Law
# =============================================================================

def gauss_law(x, y, a, g):

    G_n = 0
    N = x*y
    ind = tensor([identity(2)] * N * 2)
    vac = void(x, y)

    j = 0
    while j < y:
        i = 0
        while i < x:
            
            inner_prod = tensor(2*g*ferm_1(x, y, i, j).dag()*ferm_1(x, y, i, j) + 2*g*ferm_2(x, y, i, j).dag()*ferm_2(x, y, i, j), vac)
            E = 1/a*(Ex(x, y, i, j, a, g) - Ex(x, y, i-1, j, a, g) + Ey(x, y, i, j, a, g) - Ey(x, y, i, j-1, a, g))
            E = tensor(ind , E)
            
            G_n +=  inner_prod - E





            i += 1
        j += 1


    return G_n

def generate_projector_spin12(x, y):
#    First generate a list of all states represented by a string of 0's and 1's
    Nm = (x*y)*4
    N = x*y*2
    G_n = [0]*x*y



    all_stateReps = list(itertools.product(*[[0,1]]*Nm, *[[-1,0,1]]*Nm))
#    Then generate the corresponding values of G at site n
#               1/2*(Sz_n-1,n - Sz_n,n+1 + sz_n
    
    all_G_n = []
    
    
    for Rep in all_stateReps:
        charge = []
        j = 0
        while j < y:
            i = 0
            while i < x:
                N = (j*y + i)*2
                N1 = x*y*2
                i1 = (i-1) % x
                j1 = (j-1) % y
                N2 = (j*y + i1)*2 
                N3 = (j1*y + i)*2
                charge.extend([Rep[N] + Rep[N+1] - 1 + (Rep[N1 + N2] + Rep[N1+N3+1] - Rep[N1 + N] - Rep[N1+N+1])]) 

                
                i += 1
            j += 1
        all_G_n.append(charge)

    indices = [i for i,x in enumerate(all_G_n) if x == G_n]
    num_indices = len(indices)       

#    for i in indices:
#        print(all_stateReps[i])

#    return num_indices


#    Initialize the projector as an array to begin with
    projector_array = np.empty([2**(Nm),0])
    
#    For each "good" configuration of spins add the corresponding state vector to the projector
    for ind in indices:
        projector_array = np.append(projector_array,tensor([basis(2,x) for x in all_stateReps[ind]]).full(),axis=1)
    
#    Turn the projector into a Qobj with the correct dimensions
    projector = Qobj(projector_array)
    projector.dims = [[2]*Nm,[num_indices]]
    
    return projector



































