import matplotlib.pyplot as plt

T = [2.23, 2.24, 2.25, 2.26]
M = [0.008418, 0.009572, 0.010671, 0.007961]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 16x16, M=0)')
plt.grid(True)
plt.show()