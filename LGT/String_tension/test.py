from Phase import *
from Projector import *
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
from Operator_projection import *
from Gauss_law_advanced import *
import timeit
from Projector_advanced import *
from Utility import *
from String_projector import *

start = timeit.default_timer()


x = 3
y = 2
i = 0
j = 0
N = x*y
Nm = (x*y)*4
target = j*y + i
C = [0, 1, 0]#[0, 4/3, -1/6]#[0, 1, 0] 
a = 1
r = 1
m = -1/2
alpha = 1/137
e = np.sqrt(alpha)

#H = Hamiltonian(x, y, m, r, a, e, C)
#print(H)

save_dir = './Data/Lattice_Projection/'
Nm = (x*y)*4
G_n = [0]*x*y
filename = 'string_Spin_1_N=' + str(Nm) + '_x=' + str(x) + '_y=' +str(y)



#G = string_projector(x, y)
#save_data(G, save_dir, filename)

G = load_data(save_dir, filename)

#length = len(G)
#print(G)





"""
H_m = mass_spin(x, y, r, m, a, C)
H_E = electric_pot(x, y, a, e)
H_B = magnetic_pot(x, y, a, e)
H_Tx = momentum_x(C, x, y, a, r)    
H_Ty = momentum_y(C, x, y, a, r)    
H = a**2*(H_m + H_Tx + H_Ty)
"""

#H = Hamiltonian(x, y, m, r, a, e, C)

"""
state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]

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
result = time_evol(H , state, times, coll_ops, e_ops, args)
psi = result.states

"""
loop1 = []
loop1_img = []
loop2 = []
loop2_img = []
loop3 = []
loop3_img = []
charge = []

#charge = times
"""
l = 0
while l < num_t:
    state = psi[l]
    loop_expt1 = expect(electric(x, y, a),state)#expect(Wilson_loop(x, y, i, j),state)
#    loop_expt2 = expect(Wilson_loop2(x, y, i, j),state)
#    loop_expt3 = expect(Wilson_loop3(x, y, i, j),state)
    loop1.append(loop_expt1)
#    loop2.append(loop_expt2)
#    loop3.append(loop_expt3)
    loop1_img.append(np.imag(loop_expt1))
#    loop2_img.append(np.imag(loop_expt2))
#    loop3_img.append(np.imag(loop_expt3))
    l += 1
"""



while e < 5:
    H = Hamiltonian(x, y, m, r, a, e, C)
    state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]
    loop_expt1 = expect(Wilson_loop(x, y, i, j),state)#expect(electric_pot(x, y, a, e),state)
    loop_expt2 = expect(Wilson_line(x, y, i, j),state)
    loop_expt3 = expect(Wilson_loop3(x, y, i, j),state)
    loop1.append(loop_expt1)
    loop2.append(loop_expt2)
    loop3.append(loop_expt3)
    loop1_img.append(np.imag(loop_expt1))
    loop2_img.append(np.imag(loop_expt2))
    loop3_img.append(np.imag(loop_expt3))
    charge.append(e)
    e += 0.05


plt.figure("Time evolution")
plt.plot(charge, loop1, 'b-')
#plt.plot(charge, loop2, 'r-')
#plt.plot(charge, loop3, 'g-')
#plt.plot(charge, loop1_img, 'b--')
#plt.plot(charge, loop2_img, 'r--')
#plt.plot(charge, loop3_img, 'g--')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('Time evolution on 2x2 lattice')
plt.xlabel('time')
plt.ylabel('Number of particles')
#plt.axis([0, 5, -1, 1])
#plt.savefig('Time_evolution.svg')  
plt.show()









stop = timeit.default_timer()

print('Time: ', stop - start)  



