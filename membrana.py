"""Propiedades y transporte de agua a traves de la membrana PEM."""

import numpy as np

from config import F
from Difusividad import Difusividad
from funcion_lambda import funcion_lambda


def propiedades_membrana(actividad, T, parametros):
    """Calcula lambda, difusividad, conductividad y concentraciones de agua."""
    lambdas = funcion_lambda(actividad)
    difusividad = Difusividad(lambdas, T)
    conductividad = (0.005139 * lambdas[2] - 0.00326) * np.exp(1268.0 * (1.0 / 303.0 - 1.0 / T))
    concentraciones = (parametros["rho_membrana"] / parametros["PM_membrana"]) * lambdas[:2]
    return lambdas, difusividad, conductividad, concentraciones


def flujo_agua_membrana(corriente, lambdas, difusividad, concentraciones, parametros):
    """Flujo molar neto anodo->catodo, en mol/s."""
    nd = 2.5 * lambdas[2] / 22.0
    arrastre = nd * corriente / F
    difusion_retroceso = parametros["A_catodo"] * difusividad * (concentraciones[1] - concentraciones[0]) / parametros["e_membrana_m"]
    return arrastre - difusion_retroceso


def resistencia_membrana(conductividad_s_cm, parametros):
    """Resistencia electrica de membrana, en ohm."""
    espesor_cm = parametros["e_membrana_m"] * 100.0
    area_cm2 = parametros["A_anodo"] * 1e4
    return espesor_cm / (conductividad_s_cm * area_cm2)
