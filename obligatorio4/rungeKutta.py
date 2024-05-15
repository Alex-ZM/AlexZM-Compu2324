 
    ################################################################
    ###  PROBLEMA DE LOS TRES CUERPOS: ALGORITMO DE RUNGE-KUTTA  ###
    ################################################################

import numpy as np
import time
import os

###################################################################################################################

# DEFINICIÓN DE PARÁMETROS
G = 6.67*10**-11
M_T = 5.9736*10**24
M_L = 0.07349*10**24
d_TL = 3.844*10**8
omega = 2.6617*10**-6
R_T = 6.37816*10**6
R_L = 1.7374*10**6

m = 500
r = R_T             # Posición inicial del cohete: radio de la Tierra
phi = 90*np.pi/180  # Latitud desde la que se lanza el cohete
v = 11.2             # Módulo de la velocidad de lanzamiento de la nave (en m/s)
anguloVelocidad = 0*np.pi/180  # Ángulo de la velocidad de lanzamiento
phiPunto = 10*np.pi/180  # Velocidad angular de la nave
#pPhi = m*r**2*phiPunto

nIteraciones = 100000
h = 10

# REESCALADO DE LOS PARÁMETROS
delta = G*M_T/d_TL**3
mu = M_L/M_T
r = r/d_TL
pr = v/d_TL*np.cos(anguloVelocidad-phi)
#pPhi = pPhi/(m*d_TL**2)
pPhi = r*v*np.sin(anguloVelocidad-phi)

###################################################################################################################

def rungeKutta(h,r,phi,pr,pPhi,omega,nIteraciones,delta,mu):

    t = 0

    k_r = np.zeros(4)
    k_phi = np.zeros(4)
    k_pr = np.zeros(4)
    k_pPhi = np.zeros(4)
    
    vr = np.zeros(nIteraciones)
    vphi = np.zeros(nIteraciones)
    vpr = np.zeros(nIteraciones)
    vpPhi = np.zeros(nIteraciones)

    hamiltoniano = np.zeros(nIteraciones-1)

    posCohete = np.zeros((nIteraciones,2))
    posLuna = np.zeros((nIteraciones,2))

    vr[0] = r
    vphi[0] = phi
    vpr[0] = pr
    vpPhi[0] = pPhi

    def f_pr(pPhi,r,delta,mu,phi,omega,t):
        rPrima = np.sqrt(1+r**2-2*r*np.cos(phi-omega*t))
        return pPhi**2/r**3 - delta*(1/r**2+mu/rPrima**3*(r-np.cos(phi-omega*t)))

    def f_phi(pPhi,r):
        return pPhi/r**2
    
    def f_pPhi(delta,mu,r,phi,omega,t):
        rPrima = np.sqrt(1+r**2-2*r*np.cos(phi-omega*t))
        return -delta*mu*r/rPrima**3*np.sin(phi-omega*t)
    
    for j in range(nIteraciones-1):

        k_r[0] = h*vpr[j]
        k_phi[0] = h*f_phi(vphi[j],vr[j])
        k_pr[0] = h*f_pr(vpPhi[j],vr[j],delta,mu,vphi[j],omega,t)
        k_pPhi[0] = h*f_pPhi(delta,mu,vr[j],vphi[j],omega,t)

        k_r[1] = h*(vpr[j]+k_pr[0]/2)
        k_phi[1] = h*f_phi(vphi[j]+k_phi[0]/2,vr[j]+k_r[0]/2)
        k_pr[1] = h*f_pr(vpPhi[j]+k_phi[0]/2,vr[j]+k_r[0]/2,delta,mu,vphi[j]+k_phi[0]/2,omega,t+h/2)
        k_pPhi[1] = h*f_pPhi(delta,mu,vr[j]+k_r[0]/2,vphi[j]+k_phi[0]/2,omega,t+h/2)

        k_r[2] = h*(vpr[j]+k_pr[1]/2)
        k_phi[2] = h*f_phi(vphi[j]+k_phi[1]/2,vr[j]+k_r[1]/2)
        k_pr[2] = h*f_pr(vpPhi[j]+k_phi[1]/2,vr[j]+k_r[1]/2,delta,mu,vphi[j]+k_phi[1]/2,omega,t+h/2)
        k_pPhi[2] = h*f_pPhi(delta,mu,vr[j]+k_r[1]/2,vphi[j]+k_phi[1]/2,omega,t+h/2)

        k_r[3] = h*(vpr[j]+k_pr[2])
        k_phi[3] = h*f_phi(vphi[j]+k_phi[2],vr[j]+k_r[2])
        k_pr[3] = h*f_pr(vpPhi[j]+k_phi[2],vr[j]+k_r[2],delta,mu,vphi[j]+k_phi[2],omega,t+h)
        k_pPhi[3] = h*f_pPhi(delta,mu,vr[j]+k_r[2],vphi[j]+k_phi[2],omega,t+h)


        vr[j+1] = vr[j] + (k_r[0] + 2*k_r[1] + 2*k_r[2] + k_r[3])/6
        vphi[j+1] = vphi[j] + (k_phi[0] + 2*k_phi[1] + 2*k_phi[2] + k_phi[3])/6
        vpr[j+1] = vpr[j] + (k_pr[0] + 2*k_pr[1] + 2*k_pr[2] + k_pr[3])/6
        vpPhi[j+1] = vpPhi[j] + (k_pPhi[0] + 2*k_pPhi[1] + 2*k_pPhi[2] + k_pPhi[3])/6

        posCohete[j,0] = vr[j]*np.cos(vphi[j])
        posCohete[j,1] = vr[j]*np.sin(vphi[j])
        posLuna[j,0] = np.cos(omega*t)
        posLuna[j,1] = np.sin(omega*t)

        hamiltoniano[j] = vpr[j]**2*m*d_TL**2/2 + vpPhi[j]**2*m*d_TL**4/2 - G*m*M_T*d_TL/vr[j] - G*m*M_L*d_TL/(np.linalg.norm(np.subtract(posCohete,posLuna)))

        t = t+h

    return posCohete,posLuna,hamiltoniano

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "cohete_data.dat"     # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

posCohete,posLuna,hamiltoniano = rungeKutta(h,r,phi,pr,pPhi,omega,nIteraciones,delta,mu)

for j in range(nIteraciones):

    fichero.write(str(posCohete[j,0]) + "," + str(posCohete[j,1]) + "\n")
    fichero.write(str(posLuna[j,0]) + "," + str(posLuna[j,1]) + "\n\n")

print(hamiltoniano[0])
print(hamiltoniano[nIteraciones-2])
