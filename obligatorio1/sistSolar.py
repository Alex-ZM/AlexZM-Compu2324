# Simulación del sistema solar
import numpy as np


class Plnt:  # Clase con las propiedades de cada planeta
    def __init__(self, make, m, r, v, a, T):
        self.make = make
        self.m = m
        self.r = r
        self.v = v
        self.a = a
        self.T = T


# Definimos algunas constantes (Sistema Internacional - reescalado)
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**(-11)  # Cte de Gravitación Universal
h = 0.0003      # <---------- Paso temporal, inverso a la precisión (CAMBIAR)
nIter = 350000  # <---------- Número de iteraciones (CAMBIAR)
t = 0


def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))


# Definimos (y reescalamos) los parámetros iniciales de los planetas
sol =      Plnt(0, m=1,                       r=np.array([0,10**-15]),               v=np.array([0,0]),                 a=np.array([0,0]), T=0)
mercurio = Plnt(1, m=330.2*10**21/masaSolar,  r=np.array([57.9*10**9/UA,10**-15]),   v=np.array([0,reescalarV(47890)]), a=np.array([0,0]), T=0)
venus =    Plnt(2, m=4868.5*10**21/masaSolar, r=np.array([108.2*10**9/UA,10**-15]),  v=np.array([0,reescalarV(35030)]), a=np.array([0,0]), T=0)
tierra =   Plnt(3, m=5973.6*10**21/masaSolar, r=np.array([1,10**-15]),               v=np.array([0,reescalarV(29790)]), a=np.array([0,0]), T=0)
marte =    Plnt(4, m=641.85*10**21/masaSolar, r=np.array([227.9*10**9/UA,10**-15]),  v=np.array([0,reescalarV(24130)]), a=np.array([0,0]), T=0)
jupiter =  Plnt(5, m=1.899*10**27/masaSolar,  r=np.array([778.6*10**9/UA,10**-15]),  v=np.array([0,reescalarV(13100)]), a=np.array([0,0]), T=0)
saturno =  Plnt(6, m=0.568*10**27/masaSolar,  r=np.array([1433.5*10**9/UA,10**-15]), v=np.array([0,reescalarV(9700)]),  a=np.array([0,0]), T=0)
urano =    Plnt(7, m=0.087*10**27/masaSolar,  r=np.array([2872.5*10**9/UA,10**-15]), v=np.array([0,reescalarV(6800)]),  a=np.array([0,0]), T=0)
neptuno =  Plnt(8, m=0.102*10**27/masaSolar,  r=np.array([4495.1*10**9/UA,10**-15]), v=np.array([0,reescalarV(5400)]),  a=np.array([0,0]), T=0)
pluton =   Plnt(9, m=12.5*10**21/masaSolar,   r=np.array([5870*10**9/UA,10**-15]),   v=np.array([0,reescalarV(4700)]),  a=np.array([0,0]), T=0)

# Introducimos los planetas en un vector para acceder a ellos más fácilmente en los bucles
planeta = [sol, mercurio, venus, tierra, marte, jupiter, saturno]


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
for j in range(nIter):

    for i in range(len(planeta)):

        if planeta[i].T==0 and planeta[i].r[1]<0:
            planeta[i].T = 2*t  # Guarda el período de cada planeta

        ficheroPosiciones.write(str(planeta[i].r[0]) + ", " + str(planeta[i].r[1]) + "\n")  # Calcula entroduce las posiciones de los planetas en el fichero
    ficheroPosiciones.write("\n")  # Para separar los grupos de datos por instante temporal

    for i in range(len(planeta)):
        planeta[i].r = evR(i)  #
        planeta[i].a = evA(i)  # Avance temporal: t = t+h
        planeta[i].v = evV(i)  #
    t = t + h                  #

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPosiciones.write("# Se ha usado h = "+str(h)+"\n")
ficheroPosiciones.write("# T("+str(1)+") = "+str(planeta[1].T/planeta[3].T*365.256)+" días terrestres (vs. 87.969)\n")
ficheroPosiciones.write("# T("+str(2)+") = "+str(planeta[2].T/planeta[3].T*365.256)+" días terrestres (vs. 224.699)\n")
ficheroPosiciones.write("# T("+str(3)+") = "+str(planeta[3].T/planeta[3].T*365.256)+" días terrestres (vs. 365.256)\n")
ficheroPosiciones.write("# T("+str(4)+") = "+str(planeta[4].T/planeta[3].T*365.256)+" días terrestres (vs. 686.979)\n")
ficheroPosiciones.write("# T("+str(5)+") = "+str(planeta[5].T/planeta[3].T*365.256)+" días terrestres (vs. 4332.589)\n")
ficheroPosiciones.write("# T("+str(6)+") = "+str(planeta[6].T/planeta[3].T*365.256)+" días terrestres (vs. 10759.23)\n")

ficheroPosiciones.close()
