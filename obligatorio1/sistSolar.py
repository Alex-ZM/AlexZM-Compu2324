# Simulación del sistema solar
from planetas import Planeta
import numpy as np

# Definimos algunas constantes
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**-11  # Cte de Gravitación Universal


def reescalarT(t):  # Función para reescalar t
    return np.sqrt(G*masaSolar/UA**3)


# Definimos los parámetros iniciales de los planetas
sol =      Planeta(0, m=1,                       r=np.array(0,0),            v=(0,0), a=0)
mercurio = Planeta(1, m=330.2*10**21/masaSolar,  r=np.array(58*10**9/UA,0),  v=(0,0), a=0)
venus =    Planeta(2, m=4868.5*10**21/masaSolar, r=np.array(108*10**9/UA,0), v=(0,0), a=0)
tierra =   Planeta(3, m=5973.6*10**21/masaSolar, r=np.array(150*10**9/UA,0), v=(0,0), a=0)
venus =    Planeta(4, m=641.85*10**21/masaSolar, r=np.array(228*10**9/UA,0), v=(0,0), a=0)

# Establecemos el paso h y el tiempo t
h = 0.001  # Variar este parámetro para cambiar la precisión
t = 0


def aceleracion(i):
    return 

