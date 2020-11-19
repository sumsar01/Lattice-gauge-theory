from Phase import *
from Projector import *
from Hamiltonian_spin import *
from qutip import *
from Fields import *
from Fields_spin12 import *
from Momentum_spin import *
from Momentum import *
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

start = timeit.default_timer()


x = 2
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

"""
H_m = mass_spin(x, y, r, m, a, C_n)
H_E = electric_pot(x, y, a, g)
H_B = magnetic_pot(x, y, a, g)
H_Tx = momentum_x_spin(C_n, x, y, a, r)    
H_Ty = momentum_y_spin(C_n, x, y, a, r)    
H = a**2*(H_m + H_Tx + H_Ty)
"""
#G = generate_projector_spin12(x, y)
#G = advanced_projector(x, y)
#G = generate_projector_number(x, y, 2)
H = Hamiltonian_spin(x, y, m, r, a, e, C)
#H = Make_magnetic_pot(x, y, a, e)

#make_projection(x, y)

print(H)

#print(H)
#G = qload('Gauss_law')
#H = qload('Hamiltonian')
#H = G.dag()*H*G
#qsave(G, 'Gauss_law')
#qsave(H, 'Hamiltonian')

#state = vac_spin(x, y)
#state = H.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[1][0]
#qsave(ground, 'groundstate')
#ground = qload('groundstate')
"""
coll_ops = []   # collapse/noise operators
e_ops = []      # expectation operator
args = {'a': a, 'e': e}

#   Time evolution

max_t = 20
num_t = 100
times = np.linspace(0,max_t,num_t)
result = time_evol(H , psi, times, coll_ops, e_ops, args)
psi = result.states
"""
"""
loop1 = []
loop1_img = []
loop2 = []
loop2_img = []
loop3 = []
loop3_img = []
charge = []

a = 0.0001

while a < 20:
    H = Hamiltonian_spin(x, y, m, r, a, e, C)
    H = G.dag()*H*G
    state = H.eigenstates(sparse=True, sort='low', eigvals=5, tol=10**(-10), maxiter=100000)[1][0]
    loop_expt1 = expect(G.dag()*Wilson_loop(x, y, i, j)*G,state)
    loop_expt2 = expect(G.dag()*Wilson_loop2(x, y, i, j)*G,state)
    loop_expt3 = expect(G.dag()*Wilson_loop3(x, y, i, j)*G,state)
    loop1.append(loop_expt1)
    loop2.append(loop_expt2)
    loop3.append(loop_expt3)
    loop1_img.append(np.imag(loop_expt1))
    loop2_img.append(np.imag(loop_expt2))
    loop3_img.append(np.imag(loop_expt3))
    charge.append(a)
    a += 0.05

plt.figure("Time evolution")
plt.plot(charge, loop1, 'b-')
plt.plot(charge, loop2, 'r-')
plt.plot(charge, loop3, 'g-')
plt.plot(charge, loop1_img, 'b--')
plt.plot(charge, loop2_img, 'r--')
plt.plot(charge, loop3_img, 'g--')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
#plt.title('Time evolution on 2x2 lattice')
#plt.xlabel('time')
#plt.ylabel('Number of particles')
plt.axis([0, 20, -1, 1])
#plt.savefig('Time_evolution.svg')  
plt.show()


"""











stop = timeit.default_timer()

print('Time: ', stop - start)  



