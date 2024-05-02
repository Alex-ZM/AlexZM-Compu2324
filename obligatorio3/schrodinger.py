    ################################################
    ### RESOLUCIÓN DE LA ECUACIÓN DE SCHRÖDINGER ###
    ################################################

import numpy as np
import time
import os

###################################################################################################################

# DEFINICIÓN DE PARÁMETROS
N = 350
nCiclos = int(N/4)
lmbd = 100

skip = 3

k_0 = 2*np.pi*nCiclos/N
s = 1/(4*k_0**2)

tiempo = 300
nIter = int(np.floor(tiempo/s))

calc2N5 = int(2*N/5)
calc3N5 = int(3*N/5)
Vj = np.zeros(int(calc2N5))
VjAux = np.zeros(calc3N5-calc2N5)
for i in range(len(VjAux)):
    VjAux[i] = 0.006*lmbd*k_0**2
Vj = np.append(Vj,VjAux)
Vj = np.append(Vj,np.zeros(N-len(Vj)))

phi = np.zeros(N,dtype=np.complex64)
for j in range(len(phi)):
    phi[j] = np.exp(1j*k_0*j-8*(4*j-N)**2/N**2)
phi[0] = 0
phi[N-1] = 0

gamma = np.zeros(N,dtype=np.complex64)
alpha = np.zeros(N,dtype=np.complex64)
A_0 = np.zeros(N,dtype=np.complex64)
for j in reversed(range(1,N)):
    A_0[j] = -2+2j/s-Vj[j]
    gamma[j] = 1/(A_0[j]+alpha[j])
    alpha[j-1] = -gamma[j]

b = np.zeros(N,dtype=np.complex64)
beta = np.zeros(N,dtype=np.complex64)
ji = np.zeros(N,dtype=np.complex64)

#def calcBeta():
#    beta[N-1] = 0
#    for j in reversed(range(N-1)):
#        beta[j-1] = gamma[j]*(b[j]-beta[j])
#    
#def calcJi():
#    for j in reversed(range(N-1)):
#        ji[j] = alpha[j]

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "schrodinger_data.dat"     # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  
tIni = time.time()

for t in range(nIter):

    for j in reversed(range(N)):
        b[j] = 4j*phi[j]/s

    beta[N-1] = 0
    for j in reversed(range(N)):
        beta[j-1] = gamma[j]*(b[j]-beta[j])
    beta[0] = 0  # ¿Sobra?

    ji[0] = 0  
    for j in range(N-1):
        ji[j+1] = alpha[j]*ji[j] + beta[j]
    ji[N-1] = 0  # ¿Sobra?

    for j in range(N):
        phi[j] = ji[j]-phi[j]

    if t%skip == 0:
        norma = 0
        for j in range(N):
            fichero.write(str(j)+","+str(np.abs(phi[j]**2))+","+str(Vj[j])+"\n")
            norma += np.abs(phi[j])
        fichero.write("\n")
        print(str(norma))

tFin = time.time()
fichero.write("# Tiempo: "+str(tFin-tIni))
    