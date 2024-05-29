
    #####################################################################################
    ### SIMULACIÓN CON DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES ###
    #####################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

###################################################################################################################

datos = [[15.17967, 0.327],
         [36.65820, 0.752],
         [75.44732, 1.511],
         [129.76879, 2.574]]

temperaturas = np.array([dato[0] for dato in datos])
presiones = np.array([dato[1] for dato in datos])

pendiente, intercepto, r_value, p_value, std_err = linregress(temperaturas, presiones)
texto_leyenda = f"P = {pendiente:.4f}*T"

# Convertir temperaturas a un array de NumPy antes de la operación de multiplicación
valores_y_ajuste = pendiente * temperaturas + intercepto

plt.plot(temperaturas, valores_y_ajuste, color='red', label=f'Ajuste lineal\n{texto_leyenda}')
plt.scatter(temperaturas, presiones)
plt.xlabel('Temperatura')
plt.ylabel('Presión')
plt.title('Presión en función de la temperatura (20 partículas en red 10x10)')
plt.grid(True)
plt.legend()
plt.show()