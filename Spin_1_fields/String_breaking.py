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



e=3
H = Hamiltonian(x, y, m, r, a, e, C)


state = H.eigenstates(sparse=True, sort='low', eigvals=10, tol=10**(-10), maxiter=100000)[1][0]
state = Wilson_loop(x, y, i, j)*state


e = 7
H = Hamiltonian(x, y, m, r, a, e, C)

coll_ops = []   # collapse/noise operators
e_ops = []      # expectation operator
args = {'a': a, 'e': e}

#   Time evolution


max_t = 20
num_t = 1000
times = np.linspace(0,max_t,num_t)
result = time_evol(H , state, times, coll_ops, e_ops, args)
psi = result.states




charge = times
loop1 = []
loop2 = []
loop3 = []
loop4 = []
loop5 = []


i = 0
while i < num_t:
    state = psi[i]
    loop_expt1 = expect(num_ope(x, y),state)
    loop1.append(np.real(loop_expt1))
    loop_expt2 = expect(num_ope_site(x, y, 0, 0),state)
    loop2.append(np.real(loop_expt2))
    loop_expt3 = expect(num_ope_site(x, y, 1, 0),state)
    loop3.append(np.real(loop_expt3))    

    loop_expt4 = expect(Ex_ope(x, y, 0, 0),state)
    loop4.append(np.real(loop_expt4))    
    loop_expt5 = expect(electric(x, y, a),state)
    loop5.append(np.real(loop_expt5))    

#    loop_expt4 = expect(num_ope_site(x, y, 0, 0),state)
#    loop4.append(np.real(loop_expt4))    
#    loop_expt5 = expect(num_ope_site(x, y, 1, 1),state)
#    loop5.append(np.real(loop_expt5))    
    i += 1

#Ex_ope(x, y, i, j)
#electric(x, y, a)

plt.figure("string breaking matter fields")
plt.plot(times, loop1, 'b-')
plt.plot(charge, loop2, 'r-')
plt.plot(charge, loop3, 'g-')
#plt.plot(charge, loop4, 'k-')
#plt.plot(charge, loop5, 'm-')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('String breaking - matter fields')
plt.xlabel('time')
plt.ylabel('Number of particles')
plt.axis([0, 20, 0, 4])
plt.savefig('String_breaking_matter_fields.svg')  
plt.show()



plt.figure("string breaking gauge fields")
#plt.plot(times, loop1, 'b-')
#plt.plot(charge, loop2, 'r-')
#plt.plot(charge, loop3, 'g-')
plt.plot(charge, loop4, 'k-')
plt.plot(charge, loop5, 'm-')
#site1 = mpatches.Patch(color='blue', label='x = 1, y = 1')
#plt.legend(handles=[site1])
plt.title('String breaking - gauge fields')
plt.xlabel('time')
plt.ylabel('Electric fields')
plt.axis([0, 20, 0, 4])
plt.savefig('String_breaking_gauge_fields.svg')  
plt.show()










stop = timeit.default_timer()

print('Time: ', stop - start)  
