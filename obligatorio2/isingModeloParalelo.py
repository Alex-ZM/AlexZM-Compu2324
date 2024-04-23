import random
import numpy as np
import os
import time
import threading


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 64    # Dimensión de la cuadrícula
T = 1.44    # Temperatura T = [0,5]
t = 2000000  # Tiempo
skip = 1000

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)


def upDown(i,N):
    if i==N-1:      #
        u = i-1    #
        d = 0    #
    elif i==0:      #
        u = N-1    # Periodicidad vertical
        d = 1    #
    else:           #
        u = i-1    #
        d = i+1  #
    global up, down
    up = u
    down = d

def leftRight(j,N):
    if j==N-1:       # 
        l = j-1   # 
        r = 0    # 
    elif j==0:       # 
        l = N-1   # Periodicidad horizontal
        r = 1    # 
    else:            # 
        l = j-1   # 
        r = j+1  #
    global left, right
    left = l
    right = r


def monteCarlo(N):
    i = random.randint(0,N-1)  # Se elije una partícula aleatoria
    j = random.randint(0,N-1)  #
    
    arribajo = threading.Thread(target=upDown,args=(i,N))
    izquierderecha = threading.Thread(target=leftRight,args=(j,N))
    arribajo.start()
    izquierderecha.start()
    arribajo.join()
    izquierderecha.join()

    deltaE = 2*s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    pE = np.exp(-deltaE/T)
    if 1<pE:     # 
        p = 1    # Evaluación de p
    else:        # 
        p = pE   # 
    n = random.random()  # Número aleatorio entre 0 y 1
    return n,p,i,j


wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "ising_data.dat"           # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  


ini = time.time()
for w in range(t):

    n,p,i,j = monteCarlo(N)

    if n<p:                 # Cambio del espín (i,j)
        s[i,j] = -s[i,j]    # 
    
    # Se guarda el estado de la red en este instante
    if w%skip==0:
        fichero.write("\n")
        for i in range(N):
            fichero.write(' '.join(map(str, s[i])) + "\n")
        

fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()
