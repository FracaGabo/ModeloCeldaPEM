"""Construccion, exportacion y visualizacion de resultados."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def crear_tabla(t, concentraciones, p_parcial, corriente, fm, actividad, lambdas, dw, conductividad, e_nernst, sobrevoltajes, voltaje, voltaje_stack):
    return pd.DataFrame({
        "Tiempo (s)": t,
        "Corriente (A)": corriente,
        "H2 (mol/m3)": concentraciones[:, 0],
        "H2O_anodo (mol/m3)": concentraciones[:, 1],
        "O2 (mol/m3)": concentraciones[:, 2],
        "H2O_catodo (mol/m3)": concentraciones[:, 3],
        "N2 (mol/m3)": concentraciones[:, 4],
        "P_parcial_H2 (bar)": p_parcial[:, 0],
        "P_parcial_H2O_anodo (bar)": p_parcial[:, 1],
        "P_parcial_O2 (bar)": p_parcial[:, 2],
        "P_parcial_H2O_catodo (bar)": p_parcial[:, 3],
        "P_parcial_N2 (bar)": p_parcial[:, 4],
        "P_total_anodo (bar)": p_parcial[:, 0] + p_parcial[:, 1],
        "P_total_catodo (bar)": p_parcial[:, 2] + p_parcial[:, 3] + p_parcial[:, 4],
        "FM (mol/s)": fm,
        "Actividad_anodo": actividad[:, 0],
        "Actividad_catodo": actividad[:, 1],
        "Lambda_anodo": lambdas[:, 0],
        "Lambda_catodo": lambdas[:, 1],
        "Lambda_promedio": lambdas[:, 2],
        "Difusividad (m2/s)": dw,
        "Conductividad_membrana (S/cm)": conductividad,
        "Voltaje_Nernst (V)": e_nernst,
        "Sobrevoltaje_Anodo (V)": sobrevoltajes[:, 0],
        "Sobrevoltaje_Catodo (V)": sobrevoltajes[:, 1],
        "Sobrevoltaje_Ohmico (V)": sobrevoltajes[:, 2],
        "Sobrevoltaje_Concentracion (V)": sobrevoltajes[:, 3],
        "Voltaje_Celda (V)": voltaje,
        "Voltaje_Stack (V)": voltaje_stack,
    })


def exportar_excel(tabla, ruta):
    ruta = Path(ruta)
    tabla.to_excel(ruta, index=False, sheet_name="Datos")
    return ruta.resolve()


def crear_graficos(t, concentraciones, e_nernst, voltaje, voltaje_stack, sobrevoltajes):
    fig0, ax0 = plt.subplots(figsize=(10, 6))
    etiquetas = ("H2", "H2O (anodo)", "O2", "H2O (catodo)", "N2")
    for indice, etiqueta in enumerate(etiquetas):
        ax0.plot(t, concentraciones[:, indice], label=etiqueta, linewidth=2)
    _formatear(ax0, "Concentracion de especies vs tiempo", "Concentracion (mol/m3)", t)

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(t, e_nernst, label="Voltaje de Nernst", linewidth=2)
    ax1.plot(t, voltaje, label="Voltaje de la celda", linewidth=2)
    ax1.plot(t, voltaje_stack, label="Voltaje del stack", linewidth=2)
    _formatear(ax1, "Voltajes", "Voltaje (V)", t)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    etiquetas_eta = ("Activacion anodo", "Activacion catodo", "Ohmico", "Concentracion")
    for indice, etiqueta in enumerate(etiquetas_eta):
        ax2.plot(t, sobrevoltajes[:, indice], label=etiqueta, linewidth=2)
    _formatear(ax2, "Sobrevoltajes", "Sobrevoltaje (V)", t)
    return fig0, fig1, fig2


def _formatear(ax, titulo, ylabel, t):
    ax.set(xlabel="Tiempo (s)", ylabel=ylabel, title=titulo, xlim=(0, max(t)))
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.figure.tight_layout()
