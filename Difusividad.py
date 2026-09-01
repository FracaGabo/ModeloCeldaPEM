import numpy as np

def Difusividad(x, T):
    """Difusividad efectiva del agua en la membrana, en m2/s."""
    L = float(x[2])
    T_operacion = T

    if not np.isfinite(L) or L < 0:
        raise ValueError("El contenido de agua lambda debe ser finito y no negativo")
    if T_operacion <= 0:
        raise ValueError("La temperatura absoluta debe ser positiva")

    if L < 2:
        D_lambda = 1e-10

    elif 2 <= L <= 3:
        D_lambda = 1e-10 * (1 + 2*(L - 2))

    elif 3 < L < 4.5:
        D_lambda = 1e-10 * (3 - 1.67*(L - 3))

    elif L >= 4.5:
        D_lambda = 1.25e-10

    Dw = D_lambda * np.exp(2416*((1/303) - (1/T_operacion)))

    return Dw
