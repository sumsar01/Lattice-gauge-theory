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


x = 5
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

#H = Hamiltonian(x, y, m, r, a, e, C)
#print(H)

loop1 = []
loop1_img = []
charge = []

while e < 3:
    H = Hamiltonian(x, y, m, r, a, e, C)
    state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]
    loop_expt1 = expect(magnetic_pot(x, y, a, e),state)#expect(Wilson_loop(x, y, i, j),state)
    loop1.append(loop_expt1)
    loop1_img.append(np.imag(loop_expt1))
    charge.append(e)
    e += 0.05

plt.figure("Time evolution")
plt.plot(charge, loop1, 'b-')
plt.plot(charge, loop1_img, 'b--')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('Time evolution on 2x2 lattice')
plt.xlabel('time')
plt.ylabel('Number of particles')
#plt.axis([0, 3, -1, 1])
#plt.savefig('Time_evolution.svg')  
plt.show()











stop = timeit.default_timer()

print('Time: ', stop - start)  



