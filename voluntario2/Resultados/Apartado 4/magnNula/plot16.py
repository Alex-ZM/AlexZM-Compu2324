import matplotlib.pyplot as plt
import numpy as np

T = [1, 2, 2.5, 3, 4, 5]
M = [-1.082, -0.982, -0.733, -0.577, -0.190, -0.241]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 16x16, M=0)')
plt.grid(True)
plt.show()