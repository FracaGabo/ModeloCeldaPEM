"""Calculo de las corrientes de alimentacion del anodo y catodo."""

from config import F, R
from presion_sat import presion_sat


def calcular_entradas(parametros, corriente):
    """Devuelve caudales y concentraciones de entrada para una corriente dada."""
    T = parametros["T_operacion"]

    p_anodo = parametros["P_anodo"]
    y_h2o_an = parametros["RH_1"] * presion_sat(T, 1) / p_anodo
    y_h2_an = 1.0 - y_h2o_an
    n_h2 = parametros["exH2"] * corriente / (2.0 * F)
    n_h2o_an = n_h2 * y_h2o_an / y_h2_an

    p_catodo = parametros["P_catodo"]
    y_h2o_ca = parametros["RH_2"] * presion_sat(T, 1) / p_catodo
    y_aire_seco = 1.0 - y_h2o_ca
    y_o2 = 0.21 * y_aire_seco
    y_n2 = 0.79 * y_aire_seco
    n_o2 = parametros["exO2"] * corriente / (4.0 * F)
    n_total_ca = n_o2 / y_o2
    n_h2o_ca = n_total_ca * y_h2o_ca
    n_n2 = n_total_ca * y_n2

    return {
        "q1": (n_h2 + n_h2o_an) * R * T / (p_anodo * 1e5),
        "q2": (n_o2 + n_h2o_ca + n_n2) * R * T / (p_catodo * 1e5),
        "CH2_1": p_anodo * y_h2_an * 1e5 / (R * T),
        "CH2O_1": p_anodo * y_h2o_an * 1e5 / (R * T),
        "CO2_2": p_catodo * y_o2 * 1e5 / (R * T),
        "CH2O_2": p_catodo * y_h2o_ca * 1e5 / (R * T),
        "CN2_2": p_catodo * y_n2 * 1e5 / (R * T),
        "I": corriente,
    }
