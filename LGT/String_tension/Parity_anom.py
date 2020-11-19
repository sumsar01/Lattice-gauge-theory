from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib.patches as mpatches


# Variables

m = 1
e = 0.5*np.sqrt(m)
E0 = np.linspace(0.1, 35, 100, endpoint=True)

# Function


J_m = -m/m*e*E0/(4*np.pi)*erf(np.sqrt(np.pi*m**2/(e*E0)))
J_an = e*E0/(4*np.pi)

J_tot = J_m+J_an

# Plotting

plt.figure("Total_parity-odd_current")

plt.plot(E0,J_tot,'r-')
plt.plot(E0,J_m,'b--')
plt.plot(E0,J_an,'k--')

J_tot_leg = mpatches.Patch(color='red', label='$J_{tot}$')
J_m_leg = mpatches.Patch(color='blue', label='$J_m$')
J_an_leg = mpatches.Patch(color='black', label='$J_{an}$')
plt.legend(handles=[J_tot_leg, J_m_leg, J_an_leg])

plt.axis([0, 35,-1,1.5])
plt.title('Total parity-odd current')
plt.xlabel('Field strength ($E_c)$')
plt.ylabel('$e|E_0|/(8\pi)$')
plt.show()


