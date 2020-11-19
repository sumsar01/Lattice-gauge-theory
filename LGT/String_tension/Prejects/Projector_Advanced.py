from qutip import *
from itertools import product
import itertools as itertools
import numpy as np 


# =============================================================================
# Function that generates the projection on to the subspace of G_n
# =============================================================================
def generate_advanced_projector(x, y):
#    Lav tom liste af de gode tilstands-repræsentationer
    Nm = (x*y)*4
    G_n = [0]*x*y
    good_stateReps = []

#    Lav tom liste af de gode tilstands-indicer
#    Indicerne er tal som fortæller hvor der er partikler og antipartikler
#    Fx. fortæller listen [1,1,2,3] os, at der er: en partikel på 1. partikel-site, 
#                                                  en antipartikel på 1. antipartikel-site
#                                                  en partikel på 2. partikel-site, 
#                                                  og en antipartikel på 3. antipartikel-site.
    good_stateIndices = []
    
#    Vi beregner de gode indicer rekursivt. Denne funktion bruges til at lave rekursionen
    def loop_rec(n,i,loop_ind,indices):
        if n >= 1:
            loop_ind += 1
            for x in range(i+1,Np-(n-1)):
                for y in range(x,Np-(n-1)):
                    indices[loop_ind] = [x,y]
                    loop_rec(n-1,y,loop_ind,indices)
        else:
            good_stateIndices.append(sum(indices[:],[]))
    
#    Her foretager vi selve rekursionen
    Np = Nm//2
    for n in range(0,Np+1):
        loop_rec(n,-1,-1,[[]]*n)
    
#    Vi skal nu lave indicerne om til repræsentationer, altså lister med 0´er og 1´er
#    som repræsenterer tilstanden af både matter sites og gauge links
#    Denne funktion bruges til den omregning
    def ind2rep(state_ind):
        P_ind = [4*i for i in state_ind[::2]]
        A_ind = [4*i+2 for i in state_ind[1::2]]
        inds  = [val for pair in zip(P_ind,A_ind) for val in pair]
        
        state_rep = np.array([0,1,1,1]*Np)
        
        for ind in P_ind:
            state_rep[ind] = 1
        for ind in A_ind:
            state_rep[ind] = 0
        
        start_ind = 0
        for n,end_ind in enumerate(inds):                
            state_rep[start_ind+1:end_ind:2] = (1-(-1)**n)/2
            start_ind = end_ind
        state_rep[start_ind+1::2] = 0
        
        return list(state_rep)
#    Her foretages selve omregningen
    Pgood_stateReps = [ind2rep(state_ind) for state_ind in good_stateIndices]
    
#    Jeg har indtil videre kun fundet tilstande, hvor den første fermion specifikt er en partikel 
#    Tilstandene, hvor den første fermion er en ANTIpartikel, kan nu findes 
#    ved at lave en transformation af nogle af de førnævnte tilstande
    Agood_stateReps = []
    for stateRep in [x for x in Pgood_stateReps if x[0] == 0]:
        CC_stateRep = np.array([0,1,1,1]*Np)
        CC_stateRep[1:-2] = (1+(-1)**np.array(stateRep[3:]))/2
        
        Agood_stateReps.append(list(CC_stateRep))
    
#    Til sidst samler vi repræsentationer, som har partikler først, og repræsentationer, som har antipartikler først
    good_stateReps = Pgood_stateReps + Agood_stateReps
        
        
#    Her gemmer jeg resultatet med "save_data", som er en funktion
#    jeg selv har defineret et andet sted. Du kan nemt finde på din egen tilsvarende
    save_dir = 'Data/Nm = ' + str(Nm) + '/'
    filename = 'good_stateRepsGn' + ''.join(str(x) for x in G_n) + boundary
    save_data(good_stateReps,save_dir,filename)
    
    return good_stateReps 