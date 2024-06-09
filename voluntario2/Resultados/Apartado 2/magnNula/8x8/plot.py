import matplotlib.pyplot as plt

T = [1, 2, 2.5, 3, 4, 5]
M = [0.749, 0.686, 0.586, 0.494, 0.402, 0.361]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 8x8, M=0)')
plt.grid(True)
plt.show()