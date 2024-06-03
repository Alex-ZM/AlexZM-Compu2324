
    ####################################################################
    ###  SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI  ###
    ####################################################################

import numpy as np
import os
import time
from numba import njit
import matplotlib.pyplot as plt

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 10      # Dimensión de la cuadrícula
T = 1       # Temperatura T = [0,5]
pmc = 40000  # Número de pasos Monte Carlo
t = pmc*N**2
skip = 10*N**2

# SELECCIÓN DE LAS CONDICIONES INICIALES
# Opciones: "magnNula", "magnAleatoria"
condIni = "magnAleatoria"  

###################################################################################################################

# CREACIÓN DE LA MATRIZ DE ESPINES s
if condIni == "magnAleatoria":
    s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
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

###################################################################################################################

# CONDICIONES DE CONTORNO - BORDES VERTICALES CON ESPÍN FIJO
@njit
def condContornoV(i,N):
    if i==N-1:      #
        up = N-2    #
        down = N-1  #
    elif i==0:      #
        up = 0      # Condiciones de contorno verticales
        down = 1    # 
    else:           # 
        up = i-1    # 
        down = i+1  # 
    return up,down

# CONDICIONES DE CONTORNO - PERIODICIDAD HORIZONTAL (CILINDRO)
@njit
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
@njit
def v_pE(T,s,i,j,u,v,left1,left2,right1,right2,up1,down2):
    deltaE = 2*s[i,j]*(s[i,right1]+s[i,left1]+s[up1,j]-s[u,right2]-s[u,left2]-s[down2,v]) # Pareja vertical
    return np.exp(-deltaE/T)

# CÁLCULO DE p - CASO DE PAREJA HORIZONTAL
@njit
def h_pE(T,s,i,j,u,v,left1,right2,up1,up2,down1,down2):
    deltaE = 2*s[i,j]*(s[up1,j]+s[down1,j]+s[i,left1]-s[up2,v]-s[down2,v]-s[u,right2]) # Pareja horizontal
    return np.exp(-deltaE/T)


# MAGNETIZACIÓN DEL SISTEMA
@njit
def magnSistema(s,N):
    magn = 0
    for i in range(N):
        for j in range(N):
            magn += s[i,j]
    return magn

# MAGNETIZACIÓN PROMEDIO DE LA MITAD SUPERIOR DEL SISTEMA
@njit
def magnArriba(s,N):
    magnUp = 0
    for i in range(int((N-1)/2)):
        for j in range(N):
            magnUp += s[i,j]
    return magnUp/(N**2/2)

# MAGNETIZACIÓN PROMEDIO DE LA MITAD INFERIOR DEL SISTEMA
@njit
def magnAbajo(s,N):
    magnDown = 0
    for i in range(int((N-1)/2),N):
        for j in range(N):
            magnDown += s[i,j]
    return magnDown/(N**2/2)


# CÁLCULO DEL CALOR ESPECÍFICO
@njit
def calcCalorEspecifico(N,T,E):
    promECuad = 0
    promE = 0
    for t in range(len(E)):
        promECuad += E[t]**2
        promE += E[t]
    promECuad = promECuad/len(E)
    promE = promE/len(E)
    return (promECuad-promE**2)/(N*T)**2

# ALGORITMO DE METROPOLIS - ISING KAWASAKI
@njit
def isingKawasaki(w,flip,s,i,j,u,v,N,T):
    if s[i,j] != s[u,v]:  
        up1,down1 = condContornoV(i,N)
        up2,down2 = condContornoV(u,N)
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
@njit
def energiaMediaPorParticula(s,N):
    E = 0
    for i in range(N):
        for j in range(N):
            up,down = condContornoV(i,N)
            left,right = condContornoH(j,N)
            E += s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    return -0.5*E/N**2

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
dirDatos = "ik_data_" + str(N) + "_" + str(pmc) + "_" + str(int(skip/N**2)) + ".dat"
dirInfo = "ik_info_" + str(N) + "_" + str(pmc) + "_" + str(int(skip/N**2)) + ".dat"
fichero = open(os.path.join(wd,dirDatos), "w")
ficheroInfo = open(os.path.join(wd,dirInfo), "w")

vMagnArriba = []
vMagnAbajo = []
vEnergiaMPP = []

# BUCLE PRINCIPAL
ini = time.time()
for w in range(t):

    # Se elijen dos partículas aleatorias vecinas [ p1=(i,j) ; p2=(u,v) ]
    i = np.random.randint(1,N-2)
    j = np.random.randint(0,N)   
    flip = np.random.choice([-1,1])
    if flip == -1:  # 
        u = i + 1   # Vecinas verticales
        v = j       # 
    else:               #
        u = i           #
        if j == N-1:    # Vecinas horizontales
            v = 0       #
        else:           #
            v = j + 1   #
    
    # Una vez elegida la pareja, se decide si intercambiar los espines o no
    w = isingKawasaki(w,flip,s,i,j,u,v,N,T)

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

# REPRESENTACIÓN GRÁFICA DE LA ENERGÍA MEDIA POR PARTÍCULA 
plt.rcParams["figure.figsize"] = (7,8)

ax6 = plt.subplot(2,1,1)
plt.plot(vMagnAbajo, label="Magn. mitad inf.")
plt.plot(vMagnArriba, label="Magn. mitad sup.")
plt.axhline(y=0, color='black', linewidth=1)
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
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones ("+f"{(t/N**2):.0f}"+" pmc)")
print("|| Skip: "+ str(skip/N**2) + " pmc\n|| Calor especifico: " + f"{(calorEspecif):.7f}")
print("|| Magnetizacion del sistema: " + str(magnetizacionSistema))
print("----> Tiempo de ejecucion: "+f"{(fin-ini):.2f}"+" s")
ficheroInfo.write("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones ("+f"{(t/N**2):.0f}"+" pmc)")
ficheroInfo.write("\n|| Skip: "+ str(skip/N**2) + " pmc\n|| Calor especifico: " + f"{(calorEspecif):.7f}")
ficheroInfo.write("\n|| Magnetizacion del sistema: " + str(magnetizacionSistema))
ficheroInfo.write("\n----> Tiempo de ejecucion: "+f"{(fin-ini):.2f}"+" s\n")
fichero.close()
