"""Potencial reversible y perdidas electroquimicas de la celda."""

import numpy as np
from scipy.optimize import root_scalar

from config import F, R


def potencial_nernst(T, p_h2, p_o2):
    """Potencial reversible PEMFC; presiones parciales en bar."""
    p_h2 = np.asarray(p_h2, dtype=float)
    p_o2 = np.asarray(p_o2, dtype=float)
    if np.any(p_h2 <= 0) or np.any(p_o2 <= 0):
        raise ValueError("Nernst requiere presiones parciales positivas")
    e0_t = 1.229 - 0.85e-3 * (T - 298.15)
    return e0_t + (R * T / (2.0 * F)) * (np.log(p_h2) + 0.5 * np.log(p_o2))


def sobrevoltaje_butler_volmer(corriente, area, j0, concentracion, concentracion_ref, alpha, T):
    """Resuelve Butler-Volmer para un electrodo y devuelve eta en V."""
    if corriente < 0 or area <= 0 or j0 <= 0 or concentracion <= 0:
        raise ValueError("Parametros no fisicos en Butler-Volmer")
    if concentracion_ref <= 0 or not 0 < alpha <= 1 or T <= 0:
        raise ValueError("Referencia, alpha y temperatura deben ser validos")
    if corriente == 0:
        return 0.0
    objetivo = corriente / area
    factor = j0 * concentracion / concentracion_ref
    beta = F / (R * T)

    def residuo(eta):
        return factor * (np.exp(alpha * beta * eta) - np.exp(-(1.0 - alpha) * beta * eta)) - objetivo

    solucion = root_scalar(residuo, bracket=(0.0, 5.0), method="brentq")
    if not solucion.converged:
        raise RuntimeError("Butler-Volmer no convergio")
    return solucion.root


def sobrevoltaje_concentracion(T, densidad_corriente, densidad_limite):
    """Perdida por transporte de masa, en V."""
    if densidad_corriente < 0:
        raise ValueError("La densidad de corriente no puede ser negativa")
    if densidad_corriente >= densidad_limite:
        raise ValueError(f"J={densidad_corriente:.6g} A/m2 debe ser menor que J_limite={densidad_limite:.6g} A/m2")
    return (R * T / (2.0 * F)) * np.log(
        densidad_limite / (densidad_limite - densidad_corriente)
    )
