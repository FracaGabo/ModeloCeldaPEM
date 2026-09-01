"""Simulacion dinamica de una celda de combustible PEM."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from alimentacion import calcular_entradas
from config import R, corriente_programada, crear_parametros
from electroquimica import potencial_nernst, sobrevoltaje_butler_volmer, sobrevoltaje_concentracion
from membrana import propiedades_membrana, resistencia_membrana
from presion_sat import presion_sat
from resolver_balances import resolver_balances
from resultados import crear_graficos, crear_tabla, exportar_excel


def ejecutar_simulacion(parametros=None, ti=0.0, tf=1800.0, dt=1.0):
    """Ejecuta el modelo y devuelve todos los resultados numericos."""
    parametros = crear_parametros() if parametros is None else parametros
    corriente_base = parametros["densidad_corriente_base"] * parametros["A_anodo"]
    entradas = calcular_entradas(parametros, corriente_base)
    C_inicial = np.array([
        entradas["CH2_1"], entradas["CH2O_1"], entradas["CO2_2"],
        entradas["CH2O_2"], entradas["CN2_2"],
    ])
    t, concentraciones, fm = resolver_balances(parametros, entradas, C_inicial, ti, tf, dt)
    T = parametros["T_operacion"]
    p_parcial = concentraciones * R * T * 1e-5
    corriente = np.array([corriente_programada(instante, parametros, corriente_base) for instante in t])
    e_nernst = potencial_nernst(T, p_parcial[:, 0], p_parcial[:, 2])

    n_puntos = len(t)
    actividad = np.zeros((n_puntos, 2))
    lambdas = np.zeros((n_puntos, 3))
    dw = np.zeros(n_puntos)
    conductividad = np.zeros(n_puntos)
    sobrevoltajes = np.zeros((n_puntos, 4))
    psat = presion_sat(T, 1)

    for i in range(n_puntos):
        actividad[i] = p_parcial[i, [1, 3]] / psat
        lambdas[i], dw[i], conductividad[i], _ = propiedades_membrana(actividad[i], T, parametros)
        sobrevoltajes[i, 0] = sobrevoltaje_butler_volmer(
            corriente[i], parametros["A_anodo"], parametros["j0_anodo"],
            concentraciones[i, 0], entradas["CH2_1"], parametros["alpha_anodo"], T,
        )
        sobrevoltajes[i, 1] = sobrevoltaje_butler_volmer(
            corriente[i], parametros["A_catodo"], parametros["j0_catodo"],
            concentraciones[i, 2], entradas["CO2_2"], parametros["alpha_catodo"], T,
        )
        sobrevoltajes[i, 2] = corriente[i] * resistencia_membrana(conductividad[i], parametros)
        densidad = corriente[i] / parametros["A_anodo"]
        sobrevoltajes[i, 3] = sobrevoltaje_concentracion(T, densidad, parametros["J_limite"])

    voltaje = e_nernst - np.sum(sobrevoltajes, axis=1)
    voltaje_stack = voltaje * parametros["n"]
    if not np.all(np.isfinite(voltaje)):
        raise RuntimeError("La simulacion produjo voltajes no finitos")
    return {
        "parametros": parametros, "t": t, "concentraciones": concentraciones,
        "FM": fm, "P_parcial": p_parcial, "corriente": corriente,
        "actividad": actividad, "lambda": lambdas, "Dw": dw,
        "conductividad": conductividad, "E_nernst": e_nernst,
        "sobrevoltajes": sobrevoltajes, "voltaje": voltaje,
        "voltaje_stack": voltaje_stack,
    }


def main():
    resultados = ejecutar_simulacion()
    tabla = crear_tabla(
        resultados["t"], resultados["concentraciones"], resultados["P_parcial"],
        resultados["corriente"], resultados["FM"], resultados["actividad"],
        resultados["lambda"], resultados["Dw"], resultados["conductividad"],
        resultados["E_nernst"], resultados["sobrevoltajes"],
        resultados["voltaje"], resultados["voltaje_stack"],
    )
    ruta = exportar_excel(tabla, Path(__file__).with_name("resultados_3.xlsx"))
    crear_graficos(
        resultados["t"], resultados["concentraciones"], resultados["E_nernst"],
        resultados["voltaje"], resultados["voltaje_stack"], resultados["sobrevoltajes"],
    )
    print(f"Datos exportados a: {ruta}")
    print(f"Filas: {len(tabla)} | Columnas: {len(tabla.columns)}")
    plt.show()


if __name__ == "__main__":
    main()
