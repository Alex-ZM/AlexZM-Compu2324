import matplotlib.pyplot as plt

T = [2.23, 2.24, 2.25, 2.26]
M = [0.020335, 0.033302, 0.012970, 0.013288]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 32x32, M=0)')
plt.grid(True)
plt.show()