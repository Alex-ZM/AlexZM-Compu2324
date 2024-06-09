import matplotlib.pyplot as plt
import numpy as np

T = [1, 2, 2.5, 3, 4, 5]
M = [-1.786, -1.767, -1.102, -1.079, -0.762, -0.545]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 32x32, M=0)')
plt.grid(True)
plt.show()