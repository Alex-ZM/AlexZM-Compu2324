import matplotlib.pyplot as plt
import numpy as np

pcNumba = [[5,   12.6724],
           [10,  13.2323],
           [15,  15.8452],
           [20,  21.5361],
           [50,  38.8200],
           [70,  40.8084],
           [100, 54.2596],
           [150, 77.218],
           [200, 102.2596],
           [250, 128.1072]]

joelNumba = [[5,   12.9082],
             [10,  14.8723],
             [15,  16.5348],
             [20,  23.0102],
             [50,  39.9124],
             [70,  42.7529],
             [100, 55.5264],
             [150, 80.7124],
             [200, 107.9810],
             [250, 137.0242]]

plt.figure(figsize=(10, 6))
plt.plot([item[0] for item in pcNumba], [item[1] for item in pcNumba], label='PC (RYZEN 5 2600)')
plt.plot([item[0] for item in joelNumba], [item[1] for item in joelNumba], label='JOEL')
plt.xlim(0, 255)
plt.xticks(np.linspace(0,255, 13))
plt.title('Tiempo de ejecución en función del número de iteraciones')
plt.xlabel('Número de iteraciones (x1000)')
plt.ylabel('Tiempo de ejecución (s)')
plt.legend()
plt.show()