import numpy as np

def Difusividad(x, T):

    L = x[2]  # Tercer valor de x
    T_operacion = T

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