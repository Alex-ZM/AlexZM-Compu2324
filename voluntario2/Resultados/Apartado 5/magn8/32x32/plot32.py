import matplotlib.pyplot as plt

T = [2.2, 2.21, 2.22, 2.23, 2.24, 2.25]
M = [0.01162, 0.012488, 0.01819, 0.01449, 0.006269, 0.00633]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 32x32, M=200)')
plt.grid(True)
plt.show()