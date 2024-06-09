import matplotlib.pyplot as plt
import numpy as np

T = [1, 2, 2.5, 3, 4, 5]
M = [0.934, 0.764, 0.529, 0.338, 0.258, 0.215]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 16x16, M=8)')
plt.grid(True)
plt.show()