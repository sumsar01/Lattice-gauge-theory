from Phase import *
from Hamiltonian_spin import *
from qutip import *
from Fields import *
from Fields_spin1 import *
from Momentum_spin import *
from Potential import *
from States import *
from Time_evolution import *
from Gauss_law12 import *
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Operators_basic import *
import timeit

start = timeit.default_timer()


x = 3
y = 2
i = 0
j = 0
N = x*y
target = j*y + i
C = [0, 1, 0]#[0, 4/3, -1/6]#[0, 1, 0] 
a = 1
r = 1
m = 1
alpha = 1/137
e = np.sqrt(alpha)

H_m = mass_spin(x, y, r, m, a, C)
H_E = electric_pot(x, y, a, e)
H_B = magnetic_pot(x, y, a, e)
H_Tx = momentum_x(C, x, y, a, r)    
H_Ty = momentum_y(C, x, y, a, r)    
H = a**2*(H_m + H_Tx + H_Ty)

state = H.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[1][0]


psi = state

m = 1
e = 3
H = Hamiltonian(x, y, m, r, a, e, C)

coll_ops = []   # collapse/noise operators
e_ops = []      # expectation operator
args = {'a': a, 'e': e}

#   Time evolution
max_t = 20
num_t = 100
times = np.linspace(0,max_t,num_t)
result = time_evol(H , psi, times, coll_ops, e_ops, args)
psi = result.states

pos_00 = []
pos_10 = []
pos_20 = []
pos_01 = []
pos_11 = []
pos_21 = []
total = []


i = 0
while i < num_t:
    states = psi[i]
    pos_00.append(expect(Wilson_line(x, y, i, j),states))#expect(num_ope_site(x, y, 0, 0),states))
    pos_10.append(expect(num_ope_site(x, y, 1, 0),states))
    pos_20.append(expect(num_ope_site(x, y, 2, 0),states))
    pos_01.append(expect(num_ope_site(x, y, 0, 1),states))
    pos_11.append(expect(num_ope_site(x, y, 1, 1),states))
    pos_21.append(expect(num_ope_site(x, y, 2, 1),states))
    total.append(expect(num_ope(x, y),states))
    i += 1

plt.figure("Time evolution")
#plt.plot(times, pos_00, '-')
#plt.plot(times, pos_10, '-')
#plt.plot(times, pos_20, '-')
#plt.plot(times, pos_01, '-')
#plt.plot(times, pos_11, '-')
plt.plot(times, pos_21, '-')
#plt.plot(times, total, '-')
site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
site2 = mpatches.Patch(color='orange', label='x = 2, y = 1')
site3 = mpatches.Patch(color='green', label='x = 1, y = 2')
site4 = mpatches.Patch(color='red', label='x = 2, y = 2')
#total = mpatches.Patch(color='red', label='total particle number')
#plt.legend(handles=[site1, site2, site3, site4])
plt.title('Time evolution on 2x2 lattice')
plt.xlabel('time')
plt.ylabel('Number of particles')
#plt.axis([0, 20, 0.02, 0.04])
plt.savefig('Vacuum_states.svg')  
plt.show()

"""
plt.figure("Particle number")
plt.plot(times, total, 'r-')
total = mpatches.Patch(color='red', label='total particle number')
plt.legend(handles=[total])
plt.title('Particle numbers on 2x2 lattice')
plt.xlabel('time')
plt.ylabel('Number of particles')
#plt.axis([0, 20, 0.1, 0.15])
plt.savefig('Vacuum_number.svg')  
plt.show()
"""