import matplotlib.pyplot as plt

T = [2, 2.5, 3, 4, 5]
M = [0.909, 0.589, 0.391, 0.307, 0.288]


# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, M, 'bo', label='Datos')  # Datos marcados con círculos azules
plt.xlabel('T')
plt.ylabel('M')
plt.title('M en función de T (red 32x32, M=200)')
plt.grid(True)
plt.show()