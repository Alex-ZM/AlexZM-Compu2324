import matplotlib.pyplot as plt
import numpy as np

T = [1, 2, 2.5, 3, 4, 5]
M = [-1.509, -1.332, -1.2509, -1.114, -0.383, -0.195]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 16x16, M=8)')
plt.grid(True)
plt.show()