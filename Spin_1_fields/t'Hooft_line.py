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
loop4 = []
loop5 = []
charge = []

phase1 = np.pi
phase2 = np.pi*3/4
phase3 = np.pi/2
phase4 = np.pi/4
phase5 = 0

while e < 3:
    H = Hamiltonian(x, y, m, r, a, e, C)
    state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]
    loop_expt1 = expect(Hooft_string(x, y, i, j, phase1),state)
    loop_expt2 = expect(Hooft_string(x, y, i, j, phase2),state)
    loop_expt3 = expect(Hooft_string(x, y, i, j, phase3),state)
    loop_expt4 = expect(Hooft_string(x, y, i, j, phase4),state)
    loop_expt5 = expect(Hooft_string(x, y, i, j, phase5),state)
    loop1.append(loop_expt1)
    loop2.append(loop_expt2)
    loop3.append(loop_expt3)
    loop4.append(loop_expt4)
    loop5.append(loop_expt5)
    loop1_img.append(np.imag(loop_expt1))
    loop2_img.append(np.imag(loop_expt2))
    loop3_img.append(np.imag(loop_expt3))
    charge.append(e)
    e += 0.01




plt.figure("Time evolution")
plt.plot(charge, loop1, 'b-')
plt.plot(charge, loop2, 'r-')
plt.plot(charge, loop3, 'g-')
plt.plot(charge, loop4, 'm-')
plt.plot(charge, loop5, 'c-')
plt.plot(charge, loop1_img, 'k--')
#plt.plot(charge, loop2_img, 'r--')
#plt.plot(charge, loop3_img, 'g--')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('Spin 1 tHooft line on 2x2 lattice')
plt.xlabel('Charge (q)')
plt.ylabel('Y(L)')
plt.axis([np.sqrt(alpha), 3, -0.1, 1.1])
plt.savefig('tHooft_line.svg')  
plt.show()


stop = timeit.default_timer()

print('Time: ', stop - start)  




