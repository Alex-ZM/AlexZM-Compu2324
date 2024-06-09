import matplotlib.pyplot as plt

T = [2, 2.5, 3, 4, 5]
M = [0.775, 0.388, 0.156, 0.114, 0.086]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 32x32, M=0)')
plt.grid(True)
plt.show()