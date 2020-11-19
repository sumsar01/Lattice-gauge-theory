from Gauss_law_advanced import *
from Storage import *
from qutip import *
from itertools import product
import itertools as itertools
import numpy as np 

###############################################################################
# Making projection
###############################################################################

def make_projection(x, y):
    
    Nm = (x*y)*4
    G_n = [0]*x*y
    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Lattice_Projection/'
    filename = 'good_stateReps_N=' + str(Nm) + '_x=' + str(x) + '_y=' +str(y)
    
    if os.path.isfile(save_dir + filename + '.p'):
        good_stateReps = load_data(save_dir, filename)
    else:
        good_stateReps = advanced_projector(x, y)
        save_data(good_stateReps, save_dir, filename)
    num_good_stateReps = len(good_stateReps)

#    Initialize the projector as an array to begin with
    projector_array = np.empty([2**(Nm),0])
    
#    For each "good" configuration of spins add the corresponding state vector to the projector
    i = 0
    starttime = time.time()
    while i < num_good_stateReps: 
        progress = round(i/num_good_stateReps*100)
        speed = round(progress/(time.time() - starttime),6)
        print('\rMaking projection, progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='')
        projector_array = np.append(projector_array,tensor([basis(2,x) for x in good_stateReps[i]]).full(),axis=1)
        i += 1
    
    
#    Turn the projector into a Qobj with the correct dimensions
    projector = Qobj(projector_array)
    projector.dims = [[2]*Nm,[num_good_stateReps]]
    
    return projector

def Transform_into_G(x, y, m, a, e, operator, name):

    save_dir = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Symmetry_transformations/'
    filename = 'Gauss_law_x=' + str(x) + '_y=' +str(y)

    if os.path.isfile(save_dir + filename + '.p'):
        G = load_data(save_dir, filename)
    else:
        G = make_projection(x, y)
        save_data(G, save_dir, filename)


    save_dir2 = '/home/rasmus/Desktop/Uni/Speciale/Program/LGT/Data/Projected_operators/'
    filename2 = str(name) + '_x=' + str(x) + '_y=' + str(y) + '_m=' + str(m) + '_a=' + str(a) + '_e=' + str(e)

    if os.path.isfile(save_dir2 + filename2 + '.p'):
        operator = load_data(save_dir2, filename2)
    else:
        operator = G.dag()*operator*G
        save_data(operator, save_dir2, filename2)

    return operator













