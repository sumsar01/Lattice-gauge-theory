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

start = timeit.default_timer()

x = 2
y = 2
i = 0
j = 0
N = x*y
Nm = (x*y)*4
target = j*y + i
C = [0, 1, 0]#[0, 4/3, -1/6]#[0, 1, 0] 
a = 1
r = 1
m = 1
alpha = 1/137
e = np.sqrt(alpha)


loop1 = []
loop1_img = []
loop2 = []
loop2_img = []
loop3 = []
loop3_img = []
charge = []





while e < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e) 
    H_Tx = momentum_x(C, x, y, a, r) 
    H_Ty = momentum_y(C, x, y, a, r)
    H = a**2*(H_m + H_Tx + H_Ty + H_E)    
    state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]
    loop_expt1 = expect(electric(x, y, a),state)
    loop_expt2 = expect(Wilson_loop(x, y, i, j).dag(),state)
    loop_expt3 = expect(Wilson_loop3(x, y, i, j),state)
    loop1.append(loop_expt1)
    loop2.append(loop_expt2)
    loop3.append(loop_expt3)
    loop1_img.append(np.imag(loop_expt1))
    loop2_img.append(np.imag(loop_expt2))
    loop3_img.append(np.imag(loop_expt3))
    charge.append(e)
    e += 0.01




plt.figure("Time evolution")
plt.plot(charge, loop1, 'b-')
#plt.plot(charge, loop2, 'r--')
#plt.plot(charge, loop3, 'g-')
plt.plot(charge, loop1_img, 'k--')
#plt.plot(charge, loop2_img, 'r--')
#plt.plot(charge, loop3_img, 'g--')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('Electric field 2x2 lattice without magnetic fields')
plt.xlabel('Charge (q)')
plt.ylabel('E-field')
plt.axis([np.sqrt(alpha), 3, -1, 5])
plt.savefig('E-field_no_magnetic.svg')  
plt.show()


stop = timeit.default_timer()

print('Time: ', stop - start)  
