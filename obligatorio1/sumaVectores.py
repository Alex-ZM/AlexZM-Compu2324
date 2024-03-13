import numpy as np


# Función para sumar dos vectores
def vecSum(vect1, vect2, dim):
    vect3 = []  # Inicializamos el vector resultado

    # Convertimos en arrays numpy para poder comprobar sus dim
    v1 = np.array(vect1)
    v2 = np.array(vect2)

    # Sumamos si sus dim son iguales
    if(v1.shape == v2.shape):  
        for i in range(dim):
            vect3.append(vect1[i]+vect2[i])

        return vect3

    # Mensaje de error si las dim no son iguales
    else:
        return "Los dos vectores no tienen la misma dimensión"


# Inicializamos los vectores 1 y 2
vect1 = []
vect2 = []

# Pedimos que se introduzca la dimensión de los vectores
dim = int(input("\nIntroduce la dimensión de los vectores: "))

# Input de los elementos de los vectores
print("\nA continuación introduce los elementos del vector 1:")
for i in range(dim):
    vect1.append(float(input("|| Elemento " + str(i+1) + ": ")))
print("A continuación introduce los elementos del vector 2:")
for i in range(dim):
    vect2.append(float(input("|| Elemento " + str(i+1) + ": ")))

# Sumamos los dos vectores y mostramos el resultado.
print("\n==> Resultado v1 + v2: ")
print(vecSum(vect1, vect2, dim))
