import matplotlib.pyplot as plt
import numpy as np


class Plnt:  # Clase con las propiedades de cada planeta
    def __init__(self, make, m):
        self.make = make
        self.m = m


G = 6.67*10**(-11)  # Constante de Gravitación Universal
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol


def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))


sol =      Plnt(0, masaSolar      )  
mercurio = Plnt(1, m=330.2*10**21 ) 
venus =    Plnt(2, m=4868.5*10**21) 
tierra =   Plnt(3, m=5973.6*10**21) 
marte =    Plnt(4, m=641.85*10**21) 
jupiter =  Plnt(5, m=1.899*10**27 ) 
saturno =  Plnt(6, m=0.568*10**27 )
urano =    Plnt(7, m=0.087*10**27 )
neptuno =  Plnt(8, m=0.102*10**27 )
pluton =   Plnt(9, m=12.5*10**21  )

planeta = [sol, mercurio, venus, tierra, marte, jupiter]
T = []
V = []


def eCinetica(m, v):
    return 0.5*m*v**2


def ePotencial(m1, m2, r):
    return -G*m1*m2/r


pos = np.loadtxt("posiciones.dat")
vel = np.loadtxt("velocidades.dat")
nDatos = len(pos)

# Deshacemos los reescalados hechos en sistSolar.py
for i in range(nDatos):
    pos[i] = pos[i]*UA
    vel[i] = vel[i]/np.sqrt(UA/(G*masaSolar))

# Escribimos los vectores T y V
for t in range(0,nDatos,len(planeta)):  # itera por el número de datos de los ficheros
    for p in range(len(planeta)):
        T.append(eCinetica(planeta[p].m, np.linalg.norm(vel[t+p])))
        ePot = 0
        for i in range(len([planeta])):
            if i!=p:
                ePot += ePotencial(planeta[p].m, planeta[i].m, np.linalg.norm(pos[t+p]))
        V.append(ePot)

# Creamos e inicializamos los vectores que usaremos para el plot
tPlot = []
vPlot = []
E = []
for p in range(len(planeta)):
    tPlot.append([])
    vPlot.append([])
    E.append([])

# Guardamos los datos en tPlot y vPlot
for i in range(0, nDatos-len(planeta), len(planeta)):
    for p in range(len(planeta)):
        tPlot[p].append(T[i+p])
        vPlot[p].append(V[i+p])
print(vPlot)

# Guardamos los datos en E
for i in range(0, len(tPlot[1])):
    for p in range(len(planeta)):
        E[p].append(tPlot[p][i]+vPlot[p][i])
print(E)
# Plot
plt.plot(tPlot[1], label="T")
plt.plot(vPlot[1], label="V")
plt.plot(E[1], label="E")
plt.legend()
plt.show()
