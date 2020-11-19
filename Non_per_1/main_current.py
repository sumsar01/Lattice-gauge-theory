from Phase import *
from Hamiltonian_spin import *
from qutip import *
from Fields import *
from Momentum import *
from Potential import *
from States import *
from Time_evolution import *
from Current_x import *
from Current_y import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ================================================
# Main
# ================================================
def main():
    
    coll_ops = []   # collapse/noise operators
    e_ops = []      # expectation operator
  
    x = 2   # Number of x sites
    y = 2   # Number of y sites
    m = 0.5   # Mass
    r = 1   # Wilson term

    a = 0.25 # Lattice spacing
    e = 0.1*np.sqrt(m)   # Charge
    E_c = m**2/e    #critical field strength
    E = 1*E_c   # E-field strength    
    
#   Build Hamiltonian
    H = Hamiltonian(x, y, m, r, a)
    #H_spin = Hamiltonian_spin(x, y, m, r, a)
  
#   Build H for time evolution
    
    H = build_H(H)
    #H = build_H_uniformU(H)
   
    args = {'E': E, 'a': a, 'e': e}
   
#   Time evolution
    max_t = 50
    num_t = 100
    times = np.linspace(0,max_t,num_t)
    result = time_evol(H , vac(x, y), times, coll_ops, e_ops, args)
    states = result.states
    
    
    
    J_v = []
    J_a = []
    J_tot = []
    i = 0
    while i < len(times):
        
        J_v.append(j_vector_y(x, y, 1, 1, r, times[i], args, states[i]))
        J_a.append(j_axial_y(x, y, 1, 1, r, times[i], args, states[i]))
        J_tot.append(J_v[i] + J_a[i])    
        
        i += 1
    
    plt.figure("Current_y")
    plt.plot(times, np.absolute(J_v), 'r-')
    plt.plot(times, np.absolute(J_a), 'k-')
    plt.plot(times, np.absolute(J_tot), 'k--')
    vector = mpatches.Patch(color='red', label='$j_{v}^y$')
    axial = mpatches.Patch(color='black', label='$j_{a}^y$')
    plt.legend(handles=[axial, vector])
    #plt.axis([0, 200, 0, 0.75])
    plt.title('Current_y')
    plt.xlabel('time')
    plt.ylabel('J')
    plt.show()
    
    
 
    
main()

























