import matplotlib.pyplot as plt

T = [2.26, 2.27, 2.28]
M = [0.004252, 0.006626, 0.002079]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 48x48, M=0)')
plt.grid(True)
plt.show()