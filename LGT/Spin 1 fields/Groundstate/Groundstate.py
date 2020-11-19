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
e = np.sqrt(alpha)


G = generate_projector_number(x, y, 0)
#H = Hamiltonian_spin(x, y, m, r, a, e, C)
H_m = mass_spin(x, y, r, m, a, C)
H_E = electric_pot(x, y, a, e)
H_B = magnetic_pot(x, y, a, e)
H_Tx = momentum_x_spin(C, x, y, a, r)
H_Ty = momentum_y_spin(C, x, y, a, r)    

#H = G.dag()*H*G

Energy1 = []
Energy2 = []
charge = []

# B-field vs. E-field
"""
while e < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H1 = a**2*(H_m + H_Tx + H_Ty + H_B)
    H2 = a**2*(H_m + H_Tx + H_Ty + H_E)
    H1 = G.dag()*H1*G
    H2 = G.dag()*H2*G
    state1 = H1.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy1.append(state1)
    Energy2.append(state2)
    charge.append(e)
    e += 0.05

plt.figure("Energies_B_vs_E")
plt.plot(charge, Energy1, 'b-')
plt.plot(charge, Energy2, 'r-')
plt.axis([np.sqrt(alpha), 3, -40, 10])
Efield = mpatches.Patch(color='red', label='With E-field')
Bfield = mpatches.Patch(color='blue', label='With B-field')
plt.legend(handles=[Efield, Bfield])
plt.title('B-field vs. E-field')
plt.xlabel('Charge')
plt.ylabel('Energy')
plt.savefig('Energies_B_vs_E.svg')  
plt.show()
"""

# No B-field
"""
while e < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H1 = a**2*(H_m + H_Tx + H_Ty + H_E + H_B)
    H2 = a**2*(H_m + H_Tx + H_Ty + H_E)
    H1 = G.dag()*H1*G
    H2 = G.dag()*H2*G
    state1 = H1.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy1.append(state1)
    Energy2.append(state2)
    charge.append(e)
    e += 0.05

plt.figure("Energies_no_B")
plt.plot(charge, Energy1, 'b-')
plt.plot(charge, Energy2, 'r-')
plt.axis([np.sqrt(alpha), 3, -40, 10])
Efield = mpatches.Patch(color='red', label='Full')
Bfield = mpatches.Patch(color='blue', label='Without B-field')
plt.legend(handles=[Efield, Bfield])
plt.title('Energies without B-field')
plt.xlabel('Charge')
plt.ylabel('Energy')
plt.savefig('Energies_no_B.svg')  
plt.show()
"""
# No E-field
"""
while e < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H1 = a**2*(H_m + H_Tx + H_Ty + H_E + H_B)
    H2 = a**2*(H_m + H_Tx + H_Ty + H_B)
    H1 = G.dag()*H1*G
    H2 = G.dag()*H2*G
    state1 = H1.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy1.append(state1)
    Energy2.append(state2)
    charge.append(e)
    e += 0.05

plt.figure("Energies_no_E")
plt.plot(charge, Energy1, 'b-')
plt.plot(charge, Energy2, 'r-')
plt.axis([np.sqrt(alpha), 3, -40, 10])
Efield = mpatches.Patch(color='red', label='Without E-field')
Bfield = mpatches.Patch(color='blue', label='Full')
plt.legend(handles=[Efield, Bfield])
plt.title('Energies without E-field')
plt.xlabel('Charge')
plt.ylabel('Energy')
plt.savefig('Energies_no_E.svg')  
plt.show()
"""

# with different amounts of particles
"""
Energy0 = []
Energy2 = []
Energy4 = []
Energy6 = []
Energy8 = []
G2 = generate_projector_number(x, y, 2)
G4 = generate_projector_number(x, y, 4)
G6 = generate_projector_number(x, y, 6)
G8 = generate_projector_number(x, y, 8)

while e < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H = a**2*(H_m + H_Tx + H_Ty + H_E + H_B)
    H0 = G.dag()*H*G
    H2 = G2.dag()*H*G2
    H4 = G4.dag()*H*G4
    H6 = G6.dag()*H*G6
    H8 = G8.dag()*H*G8
    state0 = H0.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state4 = H4.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state6 = H6.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state8 = H8.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy0.append(state0)
    Energy2.append(state2)
    Energy4.append(state4)
    Energy6.append(state6)
    Energy8.append(state8)
    charge.append(e)
    e += 0.05

plt.figure("Energies_B_vs_E")
plt.plot(charge, Energy0, 'r-')
plt.plot(charge, Energy2, 'b-')
plt.plot(charge, Energy4, '-', color='green')
plt.plot(charge, Energy6, '-', color='orange')
plt.plot(charge, Energy8, '-', color='purple')
plt.axis([np.sqrt(alpha), 3, -40, 25])
particles0 = mpatches.Patch(color='red', label='0 particles')
particles2 = mpatches.Patch(color='blue', label='2 particles')
particles4 = mpatches.Patch(color='green', label='4 particles')
particles6 = mpatches.Patch(color='orange', label='6 particles')
particles8 = mpatches.Patch(color='purple', label='8 particles')
plt.legend(handles=[particles0, particles2, particles4, particles6, particles8])
plt.title('Energies with particles')
plt.xlabel('Charge')
plt.ylabel('Energy')
plt.savefig('Energies_with_particles.svg')  
plt.show()
"""

# with different amounts of particles no external field

m = 0.1
Energy0 = []
Energy2 = []
Energy4 = []
Energy6 = []
Energy8 = []
G2 = generate_projector_number(x, y, 2)
G4 = generate_projector_number(x, y, 4)
G6 = generate_projector_number(x, y, 6)
G8 = generate_projector_number(x, y, 8)

while m < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H = a**2*(H_m + H_Tx + H_Ty )
    H0 = G.dag()*H*G
    H2 = G2.dag()*H*G2
    H4 = G4.dag()*H*G4
    H6 = G6.dag()*H*G6
    H8 = G8.dag()*H*G8
    state0 = H0.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state4 = H4.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state6 = H6.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state8 = H8.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy0.append(state0)
    Energy2.append(state2)
    Energy4.append(state4)
    Energy6.append(state6)
    Energy8.append(state8)
    charge.append(m)
    m += 0.05

plt.figure("Energies_no_field")
plt.plot(charge, Energy0, 'r-')
plt.plot(charge, Energy2, 'b-')
plt.plot(charge, Energy4, '-', color='green')
plt.plot(charge, Energy6, '-', color='orange')
plt.plot(charge, Energy8, '-', color='purple')
plt.axis([np.sqrt(alpha), 3, -25, 25])
particles0 = mpatches.Patch(color='red', label='0 particles')
particles2 = mpatches.Patch(color='blue', label='2 particles')
particles4 = mpatches.Patch(color='green', label='4 particles')
particles6 = mpatches.Patch(color='orange', label='6 particles')
particles8 = mpatches.Patch(color='purple', label='8 particles')
plt.legend(handles=[particles0, particles2, particles4, particles6, particles8])
plt.title('Energies without external fields')
plt.xlabel('Mass')
plt.ylabel('Energy')
plt.savefig('Energies_no_external_fields.svg')  
plt.show()





# with different amounts of particles no external field varieng a
"""
a = 0.1
Energy0 = []
Energy2 = []
Energy4 = []
Energy6 = []
Energy8 = []
G2 = generate_projector_number(x, y, 2)
G4 = generate_projector_number(x, y, 4)
G6 = generate_projector_number(x, y, 6)
G8 = generate_projector_number(x, y, 8)

while a < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H = a**2*(H_m + H_Tx + H_Ty )
    H0 = G.dag()*H*G
    H2 = G2.dag()*H*G2
    H4 = G4.dag()*H*G4
    H6 = G6.dag()*H*G6
    H8 = G8.dag()*H*G8
    state0 = H0.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state4 = H4.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state6 = H6.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state8 = H8.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy0.append(state0)
    Energy2.append(state2)
    Energy4.append(state4)
    Energy6.append(state6)
    Energy8.append(state8)
    charge.append(a)
    a += 0.05

plt.figure("Energies_no_field")
plt.plot(charge, Energy0, 'r-')
plt.plot(charge, Energy2, 'b-')
plt.plot(charge, Energy4, '-', color='green')
plt.plot(charge, Energy6, '-', color='orange')
plt.plot(charge, Energy8, '-', color='purple')
plt.axis([np.sqrt(alpha), 3, -60, 60])
particles0 = mpatches.Patch(color='red', label='0 particles')
particles2 = mpatches.Patch(color='blue', label='2 particles')
particles4 = mpatches.Patch(color='green', label='4 particles')
particles6 = mpatches.Patch(color='orange', label='6 particles')
particles8 = mpatches.Patch(color='purple', label='8 particles')
plt.legend(handles=[particles0, particles2, particles4, particles6, particles8])
plt.title('Energies without external fields')
plt.xlabel('Lattice spacing (a)')
plt.ylabel('Energy')
plt.savefig('Energies_no_external_fields_a.svg')  
plt.show()
"""




# with different amounts of particles with external field varieng a
"""
a = 0.1
Energy0 = []
Energy2 = []
Energy4 = []
Energy6 = []
Energy8 = []
G2 = generate_projector_number(x, y, 2)
G4 = generate_projector_number(x, y, 4)
G6 = generate_projector_number(x, y, 6)
G8 = generate_projector_number(x, y, 8)

while a < 3:
    H_m = mass_spin(x, y, r, m, a, C)
    H_E = electric_pot(x, y, a, e)
    H_B = magnetic_pot(x, y, a, e)
    H_Tx = momentum_x_spin(C, x, y, a, r)
    H_Ty = momentum_y_spin(C, x, y, a, r) 
    H = a**2*(H_m + H_Tx + H_Ty + H_E + H_B )
    H0 = G.dag()*H*G
    H2 = G2.dag()*H*G2
    H4 = G4.dag()*H*G4
    H6 = G6.dag()*H*G6
    H8 = G8.dag()*H*G8
    state0 = H0.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state2 = H2.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state4 = H4.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state6 = H6.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    state8 = H8.eigenstates(sparse=True, sort='low', eigvals=0, tol=10**(-10), maxiter=100000)[0][0]
    Energy0.append(state0)
    Energy2.append(state2)
    Energy4.append(state4)
    Energy6.append(state6)
    Energy8.append(state8)
    charge.append(a)
    a += 0.05

plt.figure("Energies_no_field")
plt.plot(charge, Energy0, 'r-')
plt.plot(charge, Energy2, 'b-')
plt.plot(charge, Energy4, '-', color='green')
plt.plot(charge, Energy6, '-', color='orange')
plt.plot(charge, Energy8, '-', color='purple')
plt.axis([np.sqrt(alpha), 3, -200, 100])
particles0 = mpatches.Patch(color='red', label='0 particles')
particles2 = mpatches.Patch(color='blue', label='2 particles')
particles4 = mpatches.Patch(color='green', label='4 particles')
particles6 = mpatches.Patch(color='orange', label='6 particles')
particles8 = mpatches.Patch(color='purple', label='8 particles')
plt.legend(handles=[particles0, particles2, particles4, particles6, particles8])
plt.title('Energies with external fields')
plt.xlabel('Lattice spacing (a)')
plt.ylabel('Energy')
plt.savefig('Energies_with_external_fields_a.svg')  
plt.show()

"""

















stop = timeit.default_timer()
print('Time: ', stop - start)  