    ################################################
    ### RESOLUCIÓN DE LA ECUACIÓN DE SCHRÖDINGER ###
    ################################################

import numpy as np
import time
import os

###################################################################################################################

# DEFINICIÓN DE PARÁMETROS
N = 100
nCiclos = np.linspace(1,int(N/4), int(N/4),dtype=int)
lmbd = 1000

k_0 = []
for i in range(len(nCiclos)):
    k_0.append(2*np.pi*nCiclos[i]/N)

calc2N5 = int(2*N/5)
calc3N5 = int(3*N/5)
Vj = np.zeros(int(calc2N5))
Vj2 = np.zeros(calc3N5-calc2N5)
for i in range(len(Vj2)):
    Vj2[i] = lmbd*k_0[i]**2
Vj = np.append(Vj,Vj2)
Vj = np.append(Vj,np.zeros(N-len(Vj)))
print(Vj)
