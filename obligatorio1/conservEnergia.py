import matplotlib.pyplot as plt
import numpy as np


class Plnt:  # Clase con las propiedades de cada planeta
    def __init__(self, make, m, r, v):
        self.make = make
        self.m = m
        self.r = r
        self.v = v


G = 6.67*10**(-11)  # Constante de Gravitación Universal
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol


def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))


sol =      Plnt(0, m=1,                       r=np.array([0,10**-15]),               v=np.array([0,0])                 )  
mercurio = Plnt(1, m=330.2*10**21/masaSolar,  r=np.array([57.9*10**9/UA,10**-15]),   v=np.array([0,reescalarV(47890)]) ) 
venus =    Plnt(2, m=4868.5*10**21/masaSolar, r=np.array([108.2*10**9/UA,10**-15]),  v=np.array([0,reescalarV(35030)]) ) 
tierra =   Plnt(3, m=5973.6*10**21/masaSolar, r=np.array([1,10**-15]),               v=np.array([0,reescalarV(29790)]) ) 
marte =    Plnt(4, m=641.85*10**21/masaSolar, r=np.array([227.9*10**9/UA,10**-15]),  v=np.array([0,reescalarV(24130)]) ) 
jupiter =  Plnt(5, m=1.899*10**27/masaSolar,  r=np.array([778.6*10**9/UA,10**-15]),  v=np.array([0,reescalarV(13100)]) ) 
saturno =  Plnt(6, m=0.568*10**27/masaSolar,  r=np.array([1433.5*10**9/UA,10**-15]), v=np.array([0,reescalarV(9700)])  )
urano =    Plnt(7, m=0.087*10**27/masaSolar,  r=np.array([2872.5*10**9/UA,10**-15]), v=np.array([0,reescalarV(6800)])  )
neptuno =  Plnt(8, m=0.102*10**27/masaSolar,  r=np.array([4495.1*10**9/UA,10**-15]), v=np.array([0,reescalarV(5400)])  )
pluton =   Plnt(9, m=12.5*10**21/masaSolar,   r=np.array([5870*10**9/UA,10**-15]),   v=np.array([0,reescalarV(4700)])  )

planeta = [sol, mercurio, venus, tierra, marte, jupiter, saturno]


def T(m, v):
    return 0.5*m*v**2


def V(m1, m2, r):
    return -G*m1*m2/r


ficheroPosiciones = open("planets_data.dat", "r")
ficheroVelocidades = open("velocidades.dat", "r")

for i in ficheroPosiciones:
    n = i%len(planeta)  # Para elegir el planeta correcto

ficheroPosiciones.close()
ficheroVelocidades.close()
