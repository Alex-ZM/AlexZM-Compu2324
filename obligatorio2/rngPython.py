import random
import numpy as np
import matplotlib.pyplot as plt


def rngPython(nSamples, min, max):  # Generar e introducir 10 números aleatorios en un fichero .npy
    randomNumbers = []
    for i in range(nSamples):
        randomNumbers.append(random.randint(min, max))
    np.save("randomNumbersTest", np.array(randomNumbers))


def checkRandomness(min, max):  # Mostrar un histograma con los datos del fichero .npy
    lectura = np.load("randomNumbersTest.npy")
    plt.hist(lectura, max-min+1)
    plt.xlabel("Número")
    plt.ylabel("Frecuencia")
    plt.show()


# Input del número de int aleatorios que se gererarán
nSamples = int(input("\n¿Cuántos números aleatorios quieres generar en el fichero?: "))

# Input del intervalo de ints aleatorios que se generarán
min = 1
max = 0
print("\nA continuación introduce el intervalo de los números aleatorios generados:")
while min > max:
    min = int(input("|| Valor mínimo: "))
    max = int(input("|| Valor máximo: "))
rngPython(nSamples, min, max)  # Generamos el fichero con los parámetros dados

# Pregunta si se quiere mostrar el histograma
check = input("\n¿Mostrar histograma del fichero generado? (y/n): ")
if check == "y":
    checkRandomness(min, max)
