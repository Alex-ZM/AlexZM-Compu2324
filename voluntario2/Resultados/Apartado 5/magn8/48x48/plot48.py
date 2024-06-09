import matplotlib.pyplot as plt

T = [2.26, 2.27, 2.28]
M = [0.001407, 0.006310, 0.002588]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('Suscept. Magnética')
plt.title('Susceptibilidad magnética en función de T (red 48x48, M=8)')
plt.grid(True)
plt.show()