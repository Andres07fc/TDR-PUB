import numpy as np
import matplotlib.pyplot as plt
import math

def dilatacion_temporal(v, c = 3e8):
    """
    Esta funcion calcula la dilatación temporal

    args:
        v: velocidad relativa entre dos sistemas de referencia (m/s)
        c: velocidad de la luz en el vacio (m/s)

    return:
        factor de dilatación temporal
    """
    gamma = 1 / math.sqrt(1 - v**2/c**2)
    return gamma

# Rango de velocidades a analizar
velocidades = np.linspace(0, 0.99*3e8, 100) # Rango de valores de velocidad de 0 a 99% de c en 100 puntos

# Calcular dilatacion temporal para cada velocidad
factores_dilatacion_temporal = [dilatacion_temporal(v) for v in velocidades]

# Crear el gráfico
plt.plot(velocidades/3e8, factores_dilatacion_temporal)
plt.xlabel('Velocidad relativa (fracciones de c)')
plt.ylabel('Factor de dilatación temporal (años)')
plt.title('Dilatación temporal en función de la velocidad')
plt.grid(True)
plt.show()