from Projector_spin_1 import *
from Gauss_law_1 import *
from qutip import *
from itertools import product
import itertools as itertools
import numpy as np 
import timeit

start = timeit.default_timer()

x = 2
y = 2
n = 0

G = generate_projector_number(x, y, n)

print(G)
































stop = timeit.default_timer()

print('Time: ', stop - start)  
