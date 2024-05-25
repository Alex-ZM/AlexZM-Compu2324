import matplotlib.pyplot as plt
import numpy as np

tiempoNumba = [[2.8109],
               [3.1766],
               [3.6037],
               [4.2458],
               [4.9925],
               [5.9113]]

tiempoSinNumba = [[7.4761],
                  [16.4289],
                  [31.7434],
                  [46.0771],
                  [72.0057],
                  [91.0125]]

ax = plt.subplot(1,1,1)
plt.plot(np.linspace(2.0,7.0,6),tiempoNumba, label="Tiempo - Python + Numba")
plt.plot(np.linspace(2.0,7.0,6),tiempoSinNumba, label="Tiempo - Python")
plt.title("Tiempo de ejecución en función del número de planetas")
plt.xlabel("Número de planetas")
plt.ylabel("Tiempo de ejecución (s)")
plt.legend()
plt.show()