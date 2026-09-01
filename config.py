"""Configuracion y constantes compartidas del modelo PEMFC."""

from copy import deepcopy


R = 8.314462618  # J/(mol K)
F = 96485.33212  # C/mol


PARAMETROS_PREDETERMINADOS = {
    "alpha_anodo": 1.0,
    "alpha_catodo": 1.0,
    "j0_anodo": 2700.0,       # A/m2
    "j0_catodo": 1.0,         # A/m2
    "J_limite": 23000.0,      # A/m2
    "A_anodo": 2.5e-3,        # m2
    "A_catodo": 2.5e-3,       # m2
    "V_anodo": 0.005,         # m3
    "V_catodo": 0.01,         # m3
    "e_membrana_m": 50e-6,    # m (50 micrometros)
    "rho_membrana": 1980.0,   # kg/m3
    "PM_membrana": 1.1,       # kg/mol (1100 g/mol)
    "n": 110,
    "T_operacion": 348.15,    # K
    "exH2": 1.5,
    "exO2": 2.0,
    "P_anodo": 1.3,           # bar
    "P_catodo": 2.0,          # bar
    "RH_1": 0.5,
    "RH_2": 0.5,
    "densidad_corriente_base": 0.8e4,  # A/m2
    "alimentacion_sigue_corriente": True,
}


def crear_parametros(**cambios):
    """Crea una copia independiente de los parametros predeterminados."""
    parametros = deepcopy(PARAMETROS_PREDETERMINADOS)
    parametros.update(cambios)
    validar_parametros(parametros)
    return parametros


def validar_parametros(parametros):
    """Valida rangos y magnitudes indispensables antes de simular."""
    positivos = (
        "j0_anodo", "j0_catodo", "J_limite", "A_anodo", "A_catodo",
        "V_anodo", "V_catodo", "e_membrana_m", "rho_membrana",
        "PM_membrana", "T_operacion", "P_anodo", "P_catodo",
    )
    for nombre in positivos:
        if parametros[nombre] <= 0:
            raise ValueError(f"{nombre} debe ser positivo; recibido: {parametros[nombre]}")

    for nombre in ("RH_1", "RH_2"):
        if not 0 <= parametros[nombre] <= 1:
            raise ValueError(f"{nombre} debe estar entre 0 y 1")

    for nombre in ("alpha_anodo", "alpha_catodo"):
        if not 0 < parametros[nombre] <= 1:
            raise ValueError(f"{nombre} debe ser mayor que 0 y menor o igual que 1")

    if parametros["exH2"] < 1 or parametros["exO2"] < 1:
        raise ValueError("Los excesos de reactivo deben ser mayores o iguales que 1")


def corriente_programada(t, parametros, corriente_base=None):
    """Perfil de corriente total aplicado a la celda, en A."""
    if corriente_base is None:
        corriente_base = (
            parametros["densidad_corriente_base"] * parametros["A_anodo"]
        )
    if t < 800:
        return corriente_base
    if t < 1200:
        return 1.6e4 * parametros["A_catodo"]
    return 0.8e4 * parametros["A_catodo"]
