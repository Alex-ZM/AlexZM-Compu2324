# Simulación del sistema solar
import numpy as np


class Plnt:  # Clase con las propiedades de cada planeta
    def __init__(self, make, m, r, v, a):
        self.make = make
        self.m = m
        self.r = r
        self.v = v
        self.a = a


# Definimos algunas constantes (Sistema Internacional - reescalado)
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**(-11)  # Cte de Gravitación Universal
h = 0.01  # <--- Variar este parámetro para cambiar la precisión


def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))   # v/(UA*h*np.sqrt(G*masaSolar/(UA**3)))


# Definimos (y reescalamos) los parámetros iniciales de los planetas
sol =      Plnt(0, m=1,                       r=np.array([0,0]),             v=np.array([0,0]),                 a=np.array([0,0]))
mercurio = Plnt(1, m=330.2*10**21/masaSolar,  r=np.array([52*10**9/UA,0]),   v=np.array([0,reescalarV(47890)]), a=np.array([0,0]))
venus =    Plnt(2, m=4868.5*10**21/masaSolar, r=np.array([102*10**9/UA,0]),  v=np.array([0,reescalarV(35030)]), a=np.array([0,0]))
tierra =   Plnt(3, m=5973.6*10**21/masaSolar, r=np.array([1,0]),             v=np.array([0,reescalarV(29790)]), a=np.array([0,0]))
marte =    Plnt(4, m=641.85*10**21/masaSolar, r=np.array([222*10**9/UA,0]),  v=np.array([0,reescalarV(24130)]), a=np.array([0,0]))
jupiter =  Plnt(4, m=1.899*10**27/masaSolar,  r=np.array([775*10**9/UA,0]),  v=np.array([0,reescalarV(13100)]), a=np.array([0,0]))
saturno =  Plnt(4, m=0.568*10**27/masaSolar,  r=np.array([1430*10**9/UA,0]), v=np.array([0,reescalarV(9700)]),  a=np.array([0,0]))
urano =    Plnt(4, m=0.087*10**27/masaSolar,  r=np.array([2869*10**9/UA,0]), v=np.array([0,reescalarV(6800)]),  a=np.array([0,0]))
neptuno =  Plnt(4, m=0.102*10**27/masaSolar,  r=np.array([4491*10**9/UA,0]), v=np.array([0,reescalarV(5400)]),  a=np.array([0,0]))
pluton =   Plnt(4, m=12.5*10**21/masaSolar,   r=np.array([5865*10**9/UA,0]), v=np.array([0,reescalarV(4700)]),  a=np.array([0,0]))

# Introducimos los planetas en un vector para acceder a ellos más fácilmente en los bucles
planeta = [sol, mercurio, venus, tierra, marte, jupiter, saturno, urano, neptuno]


def a(i):  # Valor de la aceleración del planeta i en el instante actual
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal = aFinal - (planeta[j].m * np.subtract(planeta[i].r, planeta[j].r))/(np.linalg.norm(np.subtract(planeta[i].r, planeta[j].r)))**3
    return aFinal


def w(i):  #Valor de w del planeta i
    return (planeta[i].v + (h/2)*a(i))


def evR(i):  # Evolución temporal de la posición del planeta i
    return (planeta[i].r + h*w(i))


def evA(i):  # Evolución temporal de la aceleración del planeta i
    aFinal = np.array([0,0])
    for j in range(len(planeta)):
        if i != j:
            aFinal = aFinal - (planeta[j].m * np.subtract(evR(i), evR(j)))/(np.linalg.norm(np.subtract(evR(i), evR(j))))**3
    return aFinal


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return (w(i) + (h/2)*evA(i))


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
ficheroPosiciones = open("planets_data.dat", "w")
for j in range(2500):
    for i in range(len(planeta)):
        ficheroPosiciones.write(str(planeta[i].r[0]) + ", " + str(planeta[i].r[1]) + "\n")  # Calcula entroduce las posiciones de los planetas en el fichero
    ficheroPosiciones.write("\n")  # Para separar los grupos de datos por instante temporal
    for i in range(len(planeta)):
        planeta[i].r = evR(i)  #
        planeta[i].a = evA(i)  # Avance temporal: t = t+h
        planeta[i].v = evV(i)  #
ficheroPosiciones.close()

