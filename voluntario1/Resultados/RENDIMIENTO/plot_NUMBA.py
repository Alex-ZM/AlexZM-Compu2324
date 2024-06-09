import matplotlib.pyplot as plt
import numpy as np

conNumba = [[5, 12.6724],
            [10, 13.2323],
            [15, 15.8452],
            [20, 21.5361],
            [50, 38.8200],
            [70, 40.8084],
            [100, 54.2596],
            [150, 77.218],
            [200, 102.2596],
            [250, 128.1072]]

sinNumba = [[5, 44.0091],
            [10, 81.317],
            [15, 114.2991],
            [20, 158.45465],
            [25, 195.7529]]

plt.figure(figsize=(10, 6))
plt.plot([item[0] for item in conNumba], [item[1] for item in conNumba], label='con Numba')
plt.plot([item[0] for item in sinNumba], [item[1] for item in sinNumba], label='sin Numba')
plt.xlim(0, 255)
plt.xticks(np.linspace(0,255, 13))
plt.title('Tiempo de ejecución en función del número de iteraciones')
plt.xlabel('Número de iteraciones (x1000)')
plt.ylabel('Tiempo de ejecución (s)')
plt.legend()
plt.show()