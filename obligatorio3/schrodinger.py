    ################################################
    ### RESOLUCIÓN DE LA ECUACIÓN DE SCHRÖDINGER ###
    ################################################

import numpy as np
import time
import os

###################################################################################################################

# DEFINICIÓN DE PARÁMETROS (Parámetros test: N=200, lmbd=0.5, tiempo = 500)
N = 400
nCiclos = int(N/4)
lmbd = 0.5
tiempo = 500
skip = 3

k_0 = 2*np.pi*nCiclos/N
s = 1/(4*k_0**2)
nIter = int(np.floor(tiempo/s))

# VECTOR POTENCIAL
calc2N5 = int(2*N/5)
calc3N5 = int(3*N/5)
Vj = np.zeros(int(calc2N5))
VjAux = np.zeros(calc3N5-calc2N5)
for i in range(len(VjAux)):
    VjAux[i] = lmbd*k_0**2
Vj = np.append(Vj,VjAux)
Vj = np.append(Vj,np.zeros(N-len(Vj)))

###################################################################################################################

# CÁLCULO INICIAL - VECTOR PHI
phi = np.zeros(N,dtype=np.complex128)
for j in range(N):
    phi[j] = np.exp(1j*k_0*j)*np.exp(-8*(4*j-N)**2/N**2)
phi[0] = 0
phi[N-1] = 0

# CÁLCULO INICIAL - ALPHA
gamma = np.zeros(N,dtype=np.complex128)
alpha = np.zeros(N,dtype=np.complex128)
A_0 = np.zeros(N,dtype=np.complex128)
for j in range(N-1,0,-1):
    A_0[j] = -2+2j/s-Vj[j]
    gamma[j] = 1/(A_0[j]+alpha[j])
    alpha[j-1] = -gamma[j]

###################################################################################################################

# DEFINICIÓN - VECTORES b, beta, ji
b = np.zeros(N,dtype=np.complex128)
beta = np.zeros(N,dtype=np.complex128)
ji = np.zeros(N,dtype=np.complex128)

# APERTURA FICHERO
wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "schrodinger_data.dat"     # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

# BUCLE PRINCIPAL - CÁLCULO EVOLUCIÓN TEMPORAL DE PHI
tIni = time.time()
for t in range(nIter):

    for j in range(N):
        b[j] = 4j*phi[j]/s

    beta[N-1] = 0
    for j in range(N-1,1,-1):
        beta[j-1] = gamma[j]*(b[j]-beta[j])
    
    ji[0] = 0 
    for j in range(N-1):
        ji[j+1] = alpha[j]*ji[j] + beta[j]
    
    # EVOLUCIÓN DE PHI
    for j in range(N):
        phi[j] = ji[j]-phi[j]

    # GUARDA EN EL FICHERO CADA "skip" iteraciones
    if t%skip == 0:
        norma = 0
        for j in range(N):
            fichero.write(str(j)+","+str(np.abs(phi[j]**2))+","+str(Vj[j])+"\n")
            norma += np.abs(phi[j]**2)
        fichero.write("\n")
        print(str(norma))

tFin = time.time()
print("# Tiempo: "+str(tFin-tIni))
fichero.close()
    