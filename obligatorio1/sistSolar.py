# Simulación del sistema solar
from planetas import Planeta
import numpy as np

# Definimos algunas constantes
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**-11  # Cte de Gravitación Universal


def reescalarT(t):  # Función para reescalar t
    return np.sqrt(G*masaSolar/UA**3)


# Definimos (y normalizamos) los parámetros iniciales de los planetas
sol =      Planeta(0, m=1,                       r=np.array([0,0]),            v=np.array([0,0]), a=np.array([0,0]))
mercurio = Planeta(1, m=330.2*10**21/masaSolar,  r=np.array([58*10**9/UA,0]),  v=np.array([0,0]), a=np.array([0,0]))
venus =    Planeta(2, m=4868.5*10**21/masaSolar, r=np.array([108*10**9/UA,0]), v=np.array([0,0]), a=np.array([0,0]))
tierra =   Planeta(3, m=5973.6*10**21/masaSolar, r=np.array([150*10**9/UA,0]), v=np.array([0,0]), a=np.array([0,0]))
marte =    Planeta(4, m=641.85*10**21/masaSolar, r=np.array([228*10**9/UA,0]), v=np.array([0,0]), a=np.array([0,0]))
planeta = [sol, mercurio, venus, tierra, marte]

# Establecemos el paso h y el tiempo t
h = 0.001  # Variar este parámetro para cambiar la precisión
t = 0


def a(i):  # Valor de la aceleración del planeta i en el instante actual
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal -= (planeta[j].m * np.subtract(planeta[i].r, planeta[j].r))/(np.linalg.norm(np.subtract(planeta[i].r, planeta[j].r)))**3
    return aFinal


def w(i):  #Valor de w del planeta i
    return planeta[i].v + (h/2)*a(i)


def evR(i):  # Evolución temporal de la posición del planeta i
    return planeta[i].r + h*w(i)


def evA(i):  # Evolución temporal de la aceleración del planeta i
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal -= (planeta[j].m * np.subtract(evR(i), evR(j)))/(np.linalg.norm(np.subtract(evR(i), evR(j))))**3
    return aFinal


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return w(i) + (h/2)*evA(i)


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.


