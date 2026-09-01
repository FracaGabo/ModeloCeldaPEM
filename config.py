"""Configuracion y constantes compartidas del modelo PEMFC."""


R = 8.314462618  # J/(mol K)
F = 96485.33212  # C/mol


PARAMETROS_PREDETERMINADOS = {
    "alpha_anodo": 1.0,
    "alpha_catodo": 1.0,
    "j0_anodo": 2700.0,       # A/m2
    "j0_catodo": 0.1,         # A/m2
    "J_limite": 23000.0,      # A/m2
    "A_anodo": 2.5e-3,        # m2
    "A_catodo": 2.5e-3,       # m2
    "V_anodo": 0.005,         # m3
    "V_catodo": 0.01,         # m3
    "e_membrana_m": 50e-6,    # m (50 micrometros)
    "rho_membrana": 1980.0,   # kg/m3
    "PM_membrana": 1.1,       # kg/mol (1100 g/mol)
    "n": 330,
    "T_operacion": 348.15,    # K
    "exH2": 1.5,
    "exO2": 2.0,
    "P_anodo": 1.3,           # bar
    "P_catodo": 1.1,          # bar
    "RH_1": 0.5,
    "RH_2": 0.5,
    "densidad_corriente_base": 0.8e4,  # A/m2
    "alimentacion_sigue_corriente": True,
}


def crear_parametros(**cambios):
    """Crea los parametros de una simulacion."""
    parametros = PARAMETROS_PREDETERMINADOS.copy()
    parametros.update(cambios)
    return parametros


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
