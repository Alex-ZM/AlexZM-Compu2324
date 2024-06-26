
    ####################################################################
    ###  SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI  ###
    ####################################################################

import numpy as np
import os
import time
from numba import jit
import matplotlib.pyplot as plt

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 16      # Dimensión de la cuadrícula
T = 2.23       # Temperatura T = [0,5]
pmc = 200000  # Número de pasos Monte Carlo
t = pmc*N**2
skip = 200*N**2

# SELECCIÓN DE LAS CONDICIONES INICIALES
# Opciones: "magnNula", "magnAleatoria", "magnX"
condIni = "magnNula"
X = 8  # Valor de M en caso de "magnX"

###################################################################################################################

# CREACIÓN DE LA MATRIZ DE ESPINES s
# Magnetización aleatoria
if condIni == "magnAleatoria":
    s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
# Magnetización nula
elif condIni == "magnNula":
    s = np.ones((N,N)).astype(np.int8)
    for _ in range(int(N**2/2)-N):
        i = np.random.randint(1,N-1)
        j = np.random.randint(0,N-1)
        while s[i,j] == -1:
            i = np.random.randint(1,N-1)
            j = np.random.randint(0,N-1)
        s[i,j] = -1
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
# Magnetización no nula M=X
elif condIni == "magnX":
    s = np.ones((N,N)).astype(np.int8)
    for _ in range(int(N**2/2)-N):
        i = np.random.randint(1,N-1)
        j = np.random.randint(0,N-1)
        while s[i,j] == -1:
            i = np.random.randint(1,N-1)
            j = np.random.randint(0,N-1)
        s[i,j] = -1
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
    if X > 0:
        for _ in range(int(X/2)):
            i = np.random.randint(1,N-1)
            j = np.random.randint(0,N-1)
            while s[i,j] == 1:
                i = np.random.randint(1,N-1)
                j = np.random.randint(0,N-1)
            s[i,j] = 1
    else:
        for _ in range(int(X/2)):
            i = np.random.randint(1,N-1)
            j = np.random.randint(0,N-1)
            while s[i,j] == -1:
                i = np.random.randint(1,N-1)
                j = np.random.randint(0,N-1)
            s[i,j] = -1

###################################################################################################################

# CONDICIONES DE CONTORNO - PERIODICIDAD HORIZONTAL (CILINDRO)
@jit(nopython=True,fastmath=True)
def condContornoH(j,N):
    if j==N-1:       # 
        left = j-1   # 
        right = 0    # 
    elif j==0:       # 
        left = N-1   # Periodicidad horizontal
        right = 1    # 
    else:            # 
        left = j-1   # 
        right = j+1  # 
    return left,right


# CÁLCULO DE p - CASO DE PAREJA VERTICAL
@jit(nopython=True,fastmath=True)
def v_pE(T,s,i,j,u,v,left1,left2,right1,right2,up1,down2):
    deltaE = 2*s[i,j]*(s[i,right1]+s[i,left1]+s[up1,j]-s[u,right2]-s[u,left2]-s[down2,v]) # Pareja vertical
    return np.exp(-deltaE/T)

# CÁLCULO DE p - CASO DE PAREJA HORIZONTAL
@jit(nopython=True,fastmath=True)
def h_pE(T,s,i,j,u,v,left1,right2,up1,up2,down1,down2):
    deltaE = 2*s[i,j]*(s[up1,j]+s[down1,j]+s[i,left1]-s[up2,v]-s[down2,v]-s[u,right2]) # Pareja horizontal
    return np.exp(-deltaE/T)


# MAGNETIZACIÓN DEL SISTEMA
@jit(nopython=True,fastmath=True)
def magnSistema(s,N):
    magn = 0
    for i in range(N):
        for j in range(N):
            magn += s[i,j]
    return magn

# MAGNETIZACIÓN PROMEDIO DE LA MITAD SUPERIOR DEL SISTEMA
@jit(nopython=True,fastmath=True)
def magnArriba(s,N):
    magnUp = 0
    for i in range(int((N-1)/2)):
        for j in range(N):
            magnUp += s[i,j]
    return magnUp/(N**2/2)

# MAGNETIZACIÓN PROMEDIO DE LA MITAD INFERIOR DEL SISTEMA
@jit(nopython=True,fastmath=True)
def magnAbajo(s,N):
    magnDown = 0
    for i in range(int((N-1)/2),N):
        for j in range(N):
            magnDown += s[i,j]
    return magnDown/(N**2/2)


# CÁLCULO DEL CALOR ESPECÍFICO
@jit(nopython=True,fastmath=True)
def calcCalorEspecifico(N,T,E):
    promECuad = 0
    promE = 0
    for t in range(len(E)):
        promECuad += E[t]**2
        promE += E[t]
    promECuad = promECuad/len(E)
    promE = promE/len(E)
    return (promECuad-promE**2)/(T)**2


# ALGORITMO DE METROPOLIS - ISING KAWASAKI
@jit(nopython=True,fastmath=True)
def isingKawasaki(w,s,N,T):
    i = np.random.randint(1,N-2)
    j = np.random.randint(0,N)   
    flip = np.random.choice(np.array([-1,1]))
    if i == N-2:
        flip = 1
    if flip == -1:  # 
        u = i + 1   # Vecinas verticales
        v = j       # 
    else:               #
        u = i           #
        if j == N-1:    # Vecinas horizontales
            v = 0       #
        else:           #
            v = j + 1   #

    if s[i,j] != s[u,v]:  
        up1,down1 = i-1,i+1
        up2,down2 = u-1,u+1
        left1,right1 = condContornoH(j,N)
        left2,right2 = condContornoH(v,N)

        if flip == -1:
            pE = v_pE(T,s,i,j,u,v,left1,left2,right1,right2,up1,down2)  # Caso pareja vertical
        else:
            pE = h_pE(T,s,i,j,u,v,left1,right2,up1,up2,down1,down2)  # Caso pareja horizontal

        if pE>1:     # 
            p = 1    # Evaluación de p
        else:        # 
            p = pE   # 

        # Permutación de espines s1,s2 = s2,s1
        n = np.random.rand()  # Número aleatorio entre 0 y 1
        if n<p:  
            s[i,j], s[u,v] = s[u,v], s[i,j]

        return w  # Los dos espines eran diferentes -> Avanzar w
    else:
        return w-1  # Los dos espines eran iguales -> Mantener w


# ENERGÍA MEDIA POR PARTÍCULA EN EL INSTANTE t
@jit(nopython=True,fastmath=True)
def energiaMediaPorParticula(s,N):
    E = 0
    for i in range(N):
        for j in range(N):
            up,down = i-1,i+1
            left,right = condContornoH(j,N)
            E += s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    return -0.5*E/N**2

# MAGNETIZACIÓN PROMEDIO EN CADA REGIÓN
@jit(nopython=True, fastmath=True)
def calcMagnPromedio(vMagnArriba,vMagnAbajo):
    magnPromedioArriba = 0
    magnPromedioAbajo = 0
    for t in range(int(len(vMagnArriba)/2),len(vMagnArriba)):
        magnPromedioArriba += vMagnArriba[t]
        magnPromedioAbajo += vMagnAbajo[t]
    magnPromedioArriba = magnPromedioArriba/(len(vMagnArriba)/2)
    magnPromedioAbajo = magnPromedioAbajo/(len(vMagnArriba)/2)
    return magnPromedioArriba, magnPromedioAbajo

# SUSCEPTIBILIDAD MAGNÉTICA EN CADA REGIÓN
@jit(nopython=True, fastmath=True)
def calcSusceptMagnetica(vMagn,T):
    promMCuad = 0
    promM = 0
    for t in range(len(vMagn)):
        promMCuad += vMagn[t]**2
        promM += vMagn[t]
    promMCuad = promMCuad/len(vMagn)
    promM = promM/len(vMagn)
    return (promMCuad-promM**2)/T

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
dirDatos = "ik_data_T" + str(T) + "_" + str(N) + ".dat"
dirInfo = "ik_info_T" + str(T) + "_"+ str(N) + ".dat"
fichero = open(os.path.join(wd,dirDatos), "w")
ficheroInfo = open(os.path.join(wd,dirInfo), "w")

vMagnArriba = []
vMagnAbajo = []
vEnergiaMPP = []

# BUCLE PRINCIPAL
ini = time.time()
for w in range(t):
    
    # Una vez elegida la pareja, se decide si intercambiar los espines o no
    w = isingKawasaki(w,s,N,T)

    # Se guarda el estado de la red de espines y de otras magnitudes de interés
    if w%skip==0:
        fichero.write("\n")
        for i in range(N):
            fichero.write(' '.join(map(str, s[i])) + "\n")

        vMagnArriba.append(magnArriba(s,N))
        vMagnAbajo.append(magnAbajo(s,N))
        vEnergiaMPP.append(energiaMediaPorParticula(s,N))
fin = time.time()

calorEspecif = calcCalorEspecifico(N,T,vEnergiaMPP)
magnetizacionSistema = magnSistema(s,N)
magnPromedioArriba, magnPromedioAbajo = calcMagnPromedio(vMagnArriba,vMagnAbajo)
suscMagnArriba = calcSusceptMagnetica(vMagnArriba,T)
suscMagnAbajo = calcSusceptMagnetica(vMagnAbajo,T)

# REPRESENTACIÓN GRÁFICA DE LA ENERGÍA MEDIA POR PARTÍCULA 
plt.rcParams["figure.figsize"] = (7,8)

ax6 = plt.subplot(2,1,1)
plt.plot(vMagnAbajo, label="Magn. mitad inf.")
plt.plot(vMagnArriba, label="Magn. mitad sup.")
plt.axhline(y=0, color='black', linewidth=1)
plt.axhline(y=magnPromedioAbajo, color='r', linestyle='-', label="<M_inf>: "+f"{(magnPromedioAbajo):.3f}")
plt.axhline(y=magnPromedioArriba, color='r', linestyle='-', label="<M_sup>: "+f"{(magnPromedioArriba):.3f}")
plt.title("Magnetización por mitades en función de t (M_T=" + f"{magnetizacionSistema:.0f}"+")")
plt.xlabel("t")
plt.ylabel("Magnetización")
plt.legend()

ax8 = plt.subplot(2,1,2)
plt.plot(vEnergiaMPP)
plt.title("Energía promedio por partícula")
plt.xlabel("t")
plt.ylabel("E/partícula")
plt.legend()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.5)
plt.show()

###################################################################################################################

# INFORMACIÓN SOBRE LA EJECUCIÓN DEL PROGRAMA
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(pmc*N**2)+" iteraciones ("+f"{(pmc):.0f}"+" pmc)")
print("|| Skip: "+ str(skip/N**2) + " pmc\n|| Calor especifico: " + f"{(calorEspecif):.7f}")
print("|| Magnetización promedio de la mitad superior: " + f"{(magnPromedioArriba):.3f}")
print("|| Magnetización promedio de la mitad inferior: " + f"{(magnPromedioAbajo):.3f}")
print("|| Magnetizacion del sistema: " + str(magnetizacionSistema))
print("|| Susceptibilidad magnética de la mitad superior: " + f"{(suscMagnArriba):.6f}")
print("|| Susceptibilidad magnética de la mitad inferior: " + f"{(suscMagnAbajo):.6f}")
print("----> Tiempo de ejecucion: "+f"{(fin-ini):.2f}"+" s")
ficheroInfo.write("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(pmc*N**2)+" iteraciones ("+f"{(pmc):.0f}"+" pmc)")
ficheroInfo.write("\n|| Skip: "+ str(skip/N**2) + " pmc\n|| Calor especifico: " + f"{(calorEspecif):.7f}")
ficheroInfo.write("\n|| Magnetizacion promedio de la mitad superior: " + f"{(magnPromedioArriba):.3f}")
ficheroInfo.write("\n|| Magnetizacion promedio de la mitad inferior: " + f"{(magnPromedioAbajo):.3f}")
ficheroInfo.write("\n|| Magnetizacion del sistema: " + str(magnetizacionSistema))
ficheroInfo.write("\n|| Susceptibilidad magnetica de la mitad superior: " + f"{(suscMagnArriba):.6f}")
ficheroInfo.write("\n|| Susceptibilidad magnetica de la mitad inferior: " + f"{(suscMagnAbajo):.6f}")
ficheroInfo.write("\n----> Tiempo de ejecucion: "+f"{(fin-ini):.2f}"+" s\n")
fichero.close()
