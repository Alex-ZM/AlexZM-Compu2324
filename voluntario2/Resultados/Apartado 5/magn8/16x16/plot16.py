import matplotlib.pyplot as plt

T = [2.2, 2.21, 2.22, 2.23, 2.24]
M = [0.006025, 0.006520, 0.006267, 0.009389, 0.006357]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 16x16, M=8)')
plt.grid(True)
plt.show()