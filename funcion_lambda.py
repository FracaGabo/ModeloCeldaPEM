import numpy as np

def funcion_lambda(actividad):
    """Contenido de agua de Nafion para las actividades de anodo y catodo."""
    L = np.zeros(3)

    # ANODO
    if 0 <= actividad[0] <= 1:
        L[0] = 0.043 + 17.81*actividad[0] - 39.85*actividad[0]**2 + 36.0*actividad[0]**3

    elif 1 < actividad[0] <= 3:
        L[0] = 14 + 1.4*(actividad[0] - 1)

    elif actividad[0] > 3:
        L[0] = 16.8


    # CATODO
    if 0 <= actividad[1] <= 1:
        L[1] = 0.043 + 17.81*actividad[1] - 39.85*actividad[1]**2 + 36.0*actividad[1]**3

    elif 1 < actividad[1] <= 3:
        L[1] = 14 + 1.4*(actividad[1] - 1)

    elif actividad[1] > 3:
        L[1] = 16.8


    # PROMEDIO
    L[2] = (L[0] + L[1]) / 2

    return L
