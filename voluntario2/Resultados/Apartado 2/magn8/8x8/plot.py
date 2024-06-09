import matplotlib.pyplot as plt

T = [1, 2, 2.5, 3, 4, 5]
M = [0.996, 0.874, 0.778, 0.692, 0.559, 0.518]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 8x8, M=8)')
plt.grid(True)
plt.show()