"""Potencial reversible y perdidas electroquimicas de la celda."""

import numpy as np
from scipy.optimize import fsolve

from config import F, R


def potencial_nernst(T, p_h2, p_o2):
    """Potencial reversible PEMFC; presiones parciales en bar."""
    p_h2 = np.asarray(p_h2, dtype=float)
    p_o2 = np.asarray(p_o2, dtype=float)
    e0_t = 1.229 - 0.85e-3 * (T - 298.15)
    return e0_t + (R * T / (2.0 * F)) * (np.log(p_h2) + 0.5 * np.log(p_o2))


def sobrevoltaje_butler_volmer(corriente, area, j0, concentracion, concentracion_ref, alpha, T):
    """Resuelve Butler-Volmer para un electrodo y devuelve eta en V."""
    objetivo = corriente / area
    factor = j0 * concentracion / concentracion_ref
    beta = F / (R * T)

    def residuo(eta):
        return factor * (np.exp(alpha * beta * eta) - np.exp(-(1.0 - alpha) * beta * eta)) - objetivo

    return float(fsolve(residuo, 0.1)[0])


def sobrevoltaje_concentracion(T, densidad_corriente, densidad_limite):
    """Perdida por transporte de masa, en V."""
    return (R * T / (2.0 * F)) * np.log(
        densidad_limite / (densidad_limite - densidad_corriente)
    )
