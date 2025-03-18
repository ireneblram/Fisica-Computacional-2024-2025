import numpy as np

import matplotlib.pyplot as plt

# Crear un rango de valores para x
x = np.linspace(0, 2 * np.pi, 500)

# Calcular seno y coseno
y_sin = np.sin(x)
y_cos = np.cos(x)

# Crear la gráfica
plt.figure(figsize=(8, 6))
plt.plot(x, y_sin, label='Seno', color='blue')
plt.plot(x, y_cos, label='Coseno', color='red')

# Configurar etiquetas y título
plt.title('Funciones Seno y Coseno')
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
plt.legend()
plt.grid()

# Mostrar la gráfica
plt.show()