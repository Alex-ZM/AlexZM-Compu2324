import matplotlib.pyplot as plt
import numpy as np

T = [2, 2.5, 3, 4, 5]
M = [-1.739, -1.587, -0.995, -0.524, -0.312]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 32x32, M=200)')
plt.grid(True)
plt.show()