import numpy as np
import matplotlib.pyplot as plt

def contraccion_longitud(v, c = 3e8):
    """
    Esta función calcula la contracción de longitud

    args:
        L0: longitud propia del objeto (m)
        v: velocidad relativa entre dos sistemas de referencia (m/s)
        c: velocidad de la luz en el vacío (m/s)

    return:
        longitud contraída del objeto (m)
    """
    gamma = 1 / np.sqrt(1 - v**2/c**2)
    return gamma


# Rango de velocidades a analizar
velocidades = np.linspace(0, 0.99*3e8, 100)

# Calcular la contracción de longitud para cada velocidad
factor_contraccion = [contraccion_longitud(v) for v in velocidades]

# Crear el gráfico
plt.plot(velocidades/3e8, factor_contraccion)
plt.xlabel('Velocidad relativa (fracciones de c)')
plt.ylabel('factor de contracción ')
plt.title('Contracción de longitudes en función de la velocidad')
plt.grid(True)
plt.show()
