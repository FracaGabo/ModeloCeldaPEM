"""Compatibilidad con la interfaz residual historica de Butler-Volmer."""

import numpy as np

from config import F, R


def sobrevoltaje_activacion(x, u, a):
    """Devuelve los residuos de Butler-Volmer para anodo y catodo.

    x = [eta_anodo, eta_catodo]
    u = [C_H2, C_O2, temperatura, corriente]
    a = [j0_anodo, j0_catodo, C_H2_ref, C_O2_ref,
         alpha_anodo, alpha_catodo, area_anodo, area_catodo]
    """
    eta_anodo, eta_catodo = np.asarray(x, dtype=float)
    c_h2, c_o2, temperatura, corriente = np.asarray(u, dtype=float)
    valores = np.asarray(a, dtype=float)
    j0_an, j0_ca, c_h2_ref, c_o2_ref = valores[:4]
    alpha_an, alpha_ca, area_an, area_ca = valores[4:]
    if min(c_h2, c_o2, c_h2_ref, c_o2_ref, temperatura, area_an, area_ca) <= 0:
        return np.array([1e12, 1e12])

    beta = F / (R * temperatura)

    def corriente_bv(eta, alpha):
        positivo = np.clip(alpha * beta * eta, -700, 700)
        negativo = np.clip(-(1.0 - alpha) * beta * eta, -700, 700)
        return np.exp(positivo) - np.exp(negativo)

    residuo_an = j0_an * (c_h2 / c_h2_ref) * corriente_bv(eta_anodo, alpha_an) - corriente / area_an
    residuo_ca = j0_ca * (c_o2 / c_o2_ref) * corriente_bv(eta_catodo, alpha_ca) - corriente / area_ca
    return np.array([residuo_an, residuo_ca])
