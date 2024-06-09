import matplotlib.pyplot as plt

T = [1, 2, 2.5, 3, 4, 5]
M = [0.780, 0.754, 0.517, 0.313, 0.214, 0.184]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 16x16, M=0)')
plt.grid(True)
plt.show()