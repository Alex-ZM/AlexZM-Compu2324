 
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
phi = 30*np.pi/180  # Latitud desde la que se lanza el cohete
v = 200                  # Módulo de la velocidad de lanzamiento de la nave (en m/s)
phiPunto = 10*np.pi/180  # Velocidad angular de la nave
pPhi = m*r**2*phiPunto

nIteraciones = 100000
h = 0.01

# REESCALADO DE LOS PARÁMETROS
delta = G*M_T/d_TL**3
mu = M_L/M_T
r = r/d_TL
rPrima = np.sqrt(1+r**2-2*r*np.cos(phi))
v = v/d_TL
pr = v
pPhi = pPhi/(m*d_TL**2)
phiPunto = pPhi/r**2
prPunto = pPhi**2/r**3 - delta*(1/r**2+mu/rPrima**3*(r-np.cos(phi)))
pPhiPunto = -delta*mu*r/rPrima**3*np.sin(phi)

###################################################################################################################

def rungeKutta(h,pr,phiPunto,prPunto,pPhiPunto,nIteraciones):

    t = 0

    k_pr = np.zeros(4)
    k_phiPunto = np.zeros(4)
    k_prPunto = np.zeros(4)
    k_pPhiPunto = np.zeros(4)
    
    vpr = np.zeros(nIteraciones)
    vphiPunto = np.zeros(nIteraciones)
    vprPunto = np.zeros(nIteraciones)
    vpPhiPunto = np.zeros(nIteraciones)

    vpr[0] = pr
    vphiPunto[0] = phiPunto
    vprPunto[0] = prPunto
    vpPhiPunto[0] = pPhiPunto
    
    for j in range(nIteraciones-1):

        k_pr[0] = h*vpr[j]
        k_phiPunto[0] = h*vphiPunto[j]
        k_prPunto[0] = h*vprPunto[j]
        k_pPhiPunto[0] = h*vpPhiPunto[j]

        k_pr[1] = h*(vpr[j]+k_pr[0]/2)
        k_phiPunto[1] = h*(vphiPunto[j]+k_phiPunto[0]/2)
        k_prPunto[1] = h*(vprPunto[j]+k_prPunto[0]/2)
        k_pPhiPunto[1] =h*(vpPhiPunto[j]+k_pPhiPunto[0]/2)

        k_pr[2] = h*(vpr[j]+k_pr[1]/2)
        k_phiPunto[2] = h*(vphiPunto[j]+k_phiPunto[1]/2)
        k_prPunto[2] = h*(vprPunto[j]+k_prPunto[1]/2)
        k_pPhiPunto[2] =h*(vpPhiPunto[j]+k_pPhiPunto[1]/2)

        k_pr[3] = h*(vpr[j]+k_pr[2])
        k_phiPunto[3] = h*(vphiPunto[j]+k_phiPunto[2])
        k_prPunto[3] = h*(vprPunto[j]+k_prPunto[2])
        k_pPhiPunto[3] = h*(vpPhiPunto[j]+k_pPhiPunto[2])

        vpr[j+1] = vpr[j] + (k_pr[0] + 2*k_pr[1] + 2*k_pr[2] + k_pr[3])/6
        vphiPunto[j+1] = vphiPunto[j] + (k_phiPunto[0] + 2*k_phiPunto[1] + 2*k_phiPunto[2] + k_phiPunto[3])/6
        vprPunto[j+1] = vprPunto[j] + (k_prPunto[0] + 2*k_prPunto[1] + 2*k_prPunto[2] + k_prPunto[3])/6
        vpPhiPunto[j+1] = vpPhiPunto[j] + (k_pPhiPunto[0] + 2*k_pPhiPunto[1] + 2*k_pPhiPunto[2] + k_pPhiPunto[3])/6

        t = t+h

    return vpr,vphiPunto,vprPunto,vpPhiPunto

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "cohete_data.dat"     # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

vpr,vphiPunto,vprPunto,vpPhiPunto = rungeKutta(h,pr,phiPunto,prPunto,pPhiPunto,nIteraciones)

for j in range(nIteraciones):

    xCohete = vpr[j]*np.cos(vphiPunto[j])
    yCohete = vpr[j]*np.sin(vphiPunto[j])
    xLuna = np.cos(omega*j*h)
    yLuna = np.sin(omega*j*h)

    fichero.write(str(xCohete) + "," + str(yCohete) + "\n")
    fichero.write(str(xLuna) + "," + str(yLuna) + "\n\n")
