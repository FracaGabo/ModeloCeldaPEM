"""Integracion de los balances dinamicos de materia de la PEMFC."""

import numpy as np
from scipy.integrate import solve_ivp

from alimentacion import calcular_entradas
from config import F, R, corriente_programada
from membrana import flujo_agua_membrana, propiedades_membrana
from presion_sat import presion_sat


def resolver_balances(parametros, entradas, C_inicial, ti, tf, dt=1):
    """Resuelve balances y flujo de agua; conserva la API historica."""
    C_inicial = np.asarray(C_inicial, dtype=float)

    tiempos = np.arange(ti, tf + 0.5 * dt, dt)
    if tiempos[-1] < tf:
        tiempos = np.append(tiempos, tf)
    concentraciones = np.zeros((len(tiempos), 5))
    flujos_membrana = np.zeros(len(tiempos))
    concentraciones[0] = C_inicial
    corriente_base = entradas["I"]
    fm_anterior = 0.0

    for k, (t0, t1) in enumerate(zip(tiempos[:-1], tiempos[1:])):
        corriente = corriente_programada(t0, parametros, corriente_base)
        entradas_i = calcular_entradas(parametros, corriente)
        fm = fm_anterior
        for _ in range(100):
            nueva = _integrar_intervalo(t0, t1, concentraciones[k], fm, corriente, parametros, entradas_i)
            fm_calculado = _calcular_fm(nueva, corriente, parametros)
            error = abs(fm_calculado - fm)
            escala = max(abs(fm_calculado), abs(fm), 1e-12)
            if error <= 1e-10 + 1e-6 * escala:
                fm = fm_calculado
                break
            fm = 0.5 * fm + 0.5 * fm_calculado

        nueva = _integrar_intervalo(t0, t1, concentraciones[k], fm, corriente, parametros, entradas_i)
        concentraciones[k + 1] = nueva
        flujos_membrana[k + 1] = fm
        fm_anterior = fm

    return tiempos, concentraciones, flujos_membrana


def _integrar_intervalo(t0, t1, C0, fm, corriente, parametros, entradas):
    v_an, v_ca = parametros["V_anodo"], parametros["V_catodo"]
    q1, q2 = entradas["q1"], entradas["q2"]

    def balances(_t, C):
        c_h2, c_h2o_an, c_o2, c_h2o_ca, c_n2 = C
        suma_an = c_h2 + c_h2o_an
        suma_ca = c_o2 + c_h2o_ca + c_n2
        q3 = (q1 * (entradas["CH2_1"] + entradas["CH2O_1"]) - corriente / (2.0 * F) - fm) / suma_an
        q4 = (q2 * (entradas["CO2_2"] + entradas["CH2O_2"] + entradas["CN2_2"]) + corriente / (4.0 * F) + fm) / suma_ca
        return (
            q1 / v_an * entradas["CH2_1"] - q3 / v_an * c_h2 - corriente / (2.0 * F * v_an),
            q1 / v_an * entradas["CH2O_1"] - q3 / v_an * c_h2o_an - fm / v_an,
            q2 / v_ca * entradas["CO2_2"] - q4 / v_ca * c_o2 - corriente / (4.0 * F * v_ca),
            q2 / v_ca * entradas["CH2O_2"] - q4 / v_ca * c_h2o_ca + corriente / (2.0 * F * v_ca) + fm / v_ca,
            q2 / v_ca * entradas["CN2_2"] - q4 / v_ca * c_n2,
        )

    solucion = solve_ivp(balances, (t0, t1), C0, method="RK45", rtol=1e-7, atol=1e-9)
    return solucion.y[:, -1]


def _calcular_fm(C, corriente, parametros):
    psat_pa = presion_sat(parametros["T_operacion"], 1) * 1e5
    actividad = np.array([C[1], C[3]]) * R * parametros["T_operacion"] / psat_pa
    lambdas, dw, _conductividad, c_agua = propiedades_membrana(actividad, parametros["T_operacion"], parametros)
    return flujo_agua_membrana(corriente, lambdas, dw, c_agua, parametros)
