import matplotlib.pyplot as plt

T = [1, 2, 2.5, 3, 4, 5]
M = [-3.203, -2.547, -1.006, -0.811, -0.633, -0.461]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('E')
plt.title('E en función de T (red 8x8, M=0)')
plt.grid(True)
plt.show()