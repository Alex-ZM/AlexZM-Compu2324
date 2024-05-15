 
    ################################################################
    ###  PROBLEMA DE LOS TRES CUERPOS: ALGORITMO DE RUNGE-KUTTA  ###
    ################################################################

import numpy as np
import time
import os
from numba import jit
import matplotlib.pyplot as plt

###################################################################################################################

# DEFINICIÓN DE PARÁMETROS
G = 6.67*10**-11
M_T = 5.9736*10**24
M_L = 0.07349*10**24
d_TL = 3.844*10**8
omega = 2.6617*10**-6
R_T = 6.37816*10**6
R_L = 1.7374*10**6

m = 47000
r = R_T             # Posición inicial del cohete: radio de la Tierra
phi = 0*np.pi/180  # Latitud desde la que se lanza el cohete
v = 0.99999*np.sqrt(2*G*M_T/R_T)/d_TL            # Módulo de la velocidad de lanzamiento de la nave (en m/s)
anguloVelocidad = 15.455*np.pi/180      # Ángulo de la velocidad de lanzamiento
phiPunto = 10*np.pi/180  # Velocidad angular de la nave

nIteraciones = 5000000
h = 0.05
skip = 900
soloUltimoTercio = False

# REESCALADO DE LOS PARÁMETROS
delta = G*M_T/d_TL**3
mu = M_L/M_T
r = r/d_TL
pr = v*np.cos(anguloVelocidad-phi)
pPhi = r*v*np.sin(anguloVelocidad-phi)

###################################################################################################################
@jit(nopython=True,fastmath=True)
def rungeKutta(h,r,phi,pr,pPhi,omega,nIteraciones,skip,delta,mu,G,M_T,M_L,d_TL,m):

    t = 0

    k_r = np.zeros(4)
    k_phi = np.zeros(4)
    k_pr = np.zeros(4)
    k_pPhi = np.zeros(4)

    hamiltoniano = np.zeros(int(nIteraciones/skip)+1)

    posCohete = np.zeros((int(nIteraciones/skip)+1,2))
    posLuna = np.zeros((int(nIteraciones/skip)+1,2))
    

    def f_pr(pPhi,r,delta,mu,phi,omega,t):
        rPrima = np.sqrt(1+r**2-2*r*np.cos(phi-omega*t))
        return pPhi**2/r**3 - delta*(1/r**2+mu/rPrima**3*(r-np.cos(phi-omega*t)))

    def f_phi(pPhi,r):
        return pPhi/r**2
    
    def f_pPhi(delta,mu,r,phi,omega,t):
        rPrima = np.sqrt(1+r**2-2*r*np.cos(phi-omega*t))
        return -delta*mu*r/rPrima**3*np.sin(phi-omega*t)
    
    for j in range(nIteraciones-1):

        k_r[0] = h*pr
        k_phi[0] = h*f_phi(pPhi,r)
        k_pr[0] = h*f_pr(pPhi,r,delta,mu,phi,omega,t)
        k_pPhi[0] = h*f_pPhi(delta,mu,r,phi,omega,t)

        k_r[1] = h*(pr+k_pr[0]/2)
        k_phi[1] = h*f_phi(pPhi+k_pPhi[0]/2,r+k_r[0]/2)
        k_pr[1] = h*f_pr(pPhi+k_pPhi[0]/2,r+k_r[0]/2,delta,mu,phi+k_phi[0]/2,omega,t+h/2)
        k_pPhi[1] = h*f_pPhi(delta,mu,r+k_r[0]/2,phi+k_phi[0]/2,omega,t+h/2)

        k_r[2] = h*(pr+k_pr[1]/2)
        k_phi[2] = h*f_phi(pPhi+k_pPhi[1]/2,r+k_r[1]/2)
        k_pr[2] = h*f_pr(pPhi+k_pPhi[1]/2,r+k_r[1]/2,delta,mu,phi+k_phi[1]/2,omega,t+h/2)
        k_pPhi[2] = h*f_pPhi(delta,mu,r+k_r[1]/2,phi+k_phi[1]/2,omega,t+h/2)

        k_r[3] = h*(pr+k_pr[2])
        k_phi[3] = h*f_phi(pPhi+k_pPhi[2],r+k_r[2])
        k_pr[3] = h*f_pr(pPhi+k_pPhi[2],r+k_r[2],delta,mu,phi+k_phi[2],omega,t+h)
        k_pPhi[3] = h*f_pPhi(delta,mu,r+k_r[2],phi+k_phi[2],omega,t+h)


        r = r + (k_r[0] + 2*k_r[1] + 2*k_r[2] + k_r[3])/6
        phi = phi + (k_phi[0] + 2*k_phi[1] + 2*k_phi[2] + k_phi[3])/6
        pr = pr + (k_pr[0] + 2*k_pr[1] + 2*k_pr[2] + k_pr[3])/6
        pPhi = pPhi + (k_pPhi[0] + 2*k_pPhi[1] + 2*k_pPhi[2] + k_pPhi[3])/6


        if j%skip == 0:
            posCohete[int(j/skip),0] = r*np.cos(phi)
            posCohete[int(j/skip),1] = r*np.sin(phi)
            posLuna[int(j/skip),0] = np.cos(omega*t)
            posLuna[int(j/skip),1] = np.sin(omega*t)

            hamiltoniano[int(j/skip)] = (pr*m*d_TL)**2/(2*m) + (pPhi*m*d_TL**2)**2/(2*m*(r*d_TL)**2) - G*m*M_T/(r*d_TL) - G*m*M_L/(d_TL*np.linalg.norm(np.subtract(posCohete[int(j/skip)],posLuna[int(j/skip)]))) - omega*pPhi*m*d_TL**2

        t = t+h

    return posCohete,posLuna,hamiltoniano

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "cohete_data.dat"     # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

tIni = time.time()
posCohete,posLuna,hamiltoniano = rungeKutta(h,r,phi,pr,pPhi,omega,nIteraciones,skip,delta,mu,G,M_T,M_L,d_TL,m)
tFin = time.time()

if soloUltimoTercio:
    for j in range(int(2*nIteraciones/(3*skip)),int(nIteraciones/skip)):
        fichero.write(str(posCohete[j,0]) + "," + str(posCohete[j,1]) + "\n")
        fichero.write(str(posLuna[j,0]) + "," + str(posLuna[j,1]) + "\n\n")
else:
    for j in range(int(nIteraciones/skip)):
        fichero.write(str(posCohete[j,0]) + "," + str(posCohete[j,1]) + "\n")
        fichero.write(str(posLuna[j,0]) + "," + str(posLuna[j,1]) + "\n\n")

plt.plot(hamiltoniano)
plt.show()

print("Tiempo de ejecución: "+str(tFin-tIni))
