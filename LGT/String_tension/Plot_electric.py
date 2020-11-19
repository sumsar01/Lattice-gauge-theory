from Phase import *
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
import timeit

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
e = 1#np.sqrt(alpha)

#G = generate_projector_spin12(x, y)

H = electric_pot(x, y, a, e)



G = qload('Gauss_law')
#H = qload('Hamiltonian')

#qsave(G, 'Gauss_law')
#qsave(H, 'Hamiltonian')

#state = vac_spin(x, y)
state = H.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0]
#qsave(ground, 'groundstate')
#ground = qload('groundstate')


number = np.linspace(0,767,768)
number2 = [100, 200, 300, 400, 500, 600, 700]
energy = [1, 1, 1, 1, 1, 1, 1]
plt.figure("Eigenvalues")
plt.plot(number, state, 'r,')
plt.plot(number2, energy, 'ko')
plt.axis([0, 767, 0, 2])
plt.title('Electric field eigenvalues')
plt.xlabel('State number')
plt.ylabel('Energy')
theoretical = mpatches.Patch(color='black', label='Analytic value')
numeric = mpatches.Patch(color='red', label='Numerical value')
plt.legend(handles=[theoretical, numeric])
plt.savefig('Electric_Eigenvalues.svg')  
plt.show()