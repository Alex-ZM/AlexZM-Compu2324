import matplotlib.pyplot as plt

T = [1, 2, 2.5, 3, 4, 5]
M = [-3.141, -1.574, -1.025, -0.851, -0.603, -0.319]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 8x8, M=8)')
plt.grid(True)
plt.show()