# Simulación del sistema solar
import numpy as np


class Planeta:  # Clase con las propiedades de cada planeta
    def __init__(self, make, m, r, v, a):
        self.make = make
        self.m = m
        self.r = r
        self.v = v
        self.a = a


# Definimos algunas constantes
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**-11  # Cte de Gravitación Universal
h = 0.01  # <--- Variar este parámetro para cambiar la precisión


def reescalarV(v):  # Función para reescalar t
    return v/(UA*np.sqrt(G*masaSolar/UA**3)*h)  # C O M P R O B A R


# Definimos (y normalizamos) los parámetros iniciales de los planetas
sol =      Planeta(0, m=1,                       r=np.array([0,0]),            v=np.array([0,0]),                 a=np.array([0,0]))
mercurio = Planeta(1, m=330.2*10**21/masaSolar,  r=np.array([58*10**9/UA,0]),  v=np.array([0,reescalarV(47890)]), a=np.array([0,0]))
venus =    Planeta(2, m=4868.5*10**21/masaSolar, r=np.array([108*10**9/UA,0]), v=np.array([0,reescalarV(35030)]), a=np.array([0,0]))
tierra =   Planeta(3, m=5973.6*10**21/masaSolar, r=np.array([150*10**9/UA,0]), v=np.array([0,reescalarV(29790)]), a=np.array([0,0]))
marte =    Planeta(4, m=641.85*10**21/masaSolar, r=np.array([228*10**9/UA,0]), v=np.array([0,reescalarV(24130)]), a=np.array([0,0]))

# Introducimos los planetas en un vector para acceder a ellos más fácilmente en los bucles
planeta = [sol, mercurio, venus, tierra, marte]  # No confundir con los objetos de la clase planetas (P mayúscula)



def a(i):  # Valor de la aceleración del planeta i en el instante actual
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal = aFinal - (planeta[j].m * np.subtract(planeta[i].r, planeta[j].r))/(np.linalg.norm(np.subtract(planeta[i].r, planeta[j].r)))**3
    return aFinal


def w(i):  #Valor de w del planeta i
    return planeta[i].v + (h/2)*a(i)


def evR(i):  # Evolución temporal de la posición del planeta i
    return planeta[i].r + h*w(i)


def evA(i):  # Evolución temporal de la aceleración del planeta i
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal = aFinal - (planeta[j].m * np.subtract(evR(i), evR(j)))/(np.linalg.norm(np.subtract(evR(i), evR(j))))**3
    return aFinal


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return w(i) + (h/2)*evA(i)


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
ficheroPosiciones = open("planets_data.dat", "w")
for j in range(1000):
    for i in range(len(planeta)):
        ficheroPosiciones.write(str(planeta[i].r[0]) + ", " + str(planeta[i].r[1]) + "\n")  # Introduce las posiciones en el fichero
    ficheroPosiciones.write("\n")  # Para separar los grupos de datos por instante temporal
    for i in range(len(planeta)):
        planeta[i].r = evR(i)  #
        planeta[i].a = evA(i)  # Avance temporal: t = t+h
        planeta[i].v = evV(i)  #
ficheroPosiciones.close()

