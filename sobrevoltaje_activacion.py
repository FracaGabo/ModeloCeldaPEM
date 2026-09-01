import numpy as np

def sobrevoltaje_activacion(x, u, a):
    """
    Calcula los sobrevoltajes de activación del ánodo y cátodo.
    
    Parámetros:
    -----------
    x : array-like
        [n_anodo, n_catodo] - sobrevoltajes a resolver
    u : array-like
        [CH2, CO2, T_operacion, I] - variables de decisión
    a : array-like
        [j0_anodo, j0_catodo, CH2_1, CO2_2, alpha_anodo, alpha_catodo, A_anodo, A_catodo]
    
    Retorna:
    --------
    V_act : array
        [ecuacion_anodo, ecuacion_catodo] - sistema de ecuaciones a resolver
    """
    
    n = 1
    
    # CONSTANTES
    R = 8.314  # J/mol K
    F = 96485  # C/mol
    
    # VARIABLES DE x
    n_anodo = x[0]
    n_catodo = x[1]
    
    # VARIABLES DE u
    CH2 = u[0]
    CO2 = u[1]
    T_operacion = u[2]
    I = u[3]
    
    # VARIABLES DE a
    j0_anodo = a[0]
    j0_catodo = a[1]
    CH2_1 = a[2]
    CO2_2 = a[3]
    alpha_anodo = a[4]
    alpha_catodo = a[5]
    A_anodo = a[6]
    A_catodo = a[7]
    
    # DENSIDADES DE CORRIENTE
    J_anodo = I / A_anodo
    J_catodo = I / A_catodo
    
    # FACTORES
    f_anodo = (alpha_anodo * F) / (R * T_operacion)
    f_catodo = (alpha_catodo * F) / (R * T_operacion)
    
    # LIMITAR EXPONENTES PARA EVITAR OVERFLOW
    MAX_EXP = 700  # Límite seguro para np.exp()
    
    # ECUACIÓN 1: ÁNODO
    exp_arg_1_pos = f_anodo * n_anodo
    exp_arg_1_neg = -f_catodo * n_anodo
    
    exp_arg_1_pos = np.clip(exp_arg_1_pos, -MAX_EXP, MAX_EXP)
    exp_arg_1_neg = np.clip(exp_arg_1_neg, -MAX_EXP, MAX_EXP)
    
    V_act_1 = j0_anodo * ((CH2 / CH2_1) ** n) * (np.exp(exp_arg_1_pos) - np.exp(exp_arg_1_neg)) - J_anodo
    
    # ECUACIÓN 2: CÁTODO
    exp_arg_2_pos = f_catodo * n_catodo
    exp_arg_2_neg = -f_anodo * n_catodo
    
    exp_arg_2_pos = np.clip(exp_arg_2_pos, -MAX_EXP, MAX_EXP)
    exp_arg_2_neg = np.clip(exp_arg_2_neg, -MAX_EXP, MAX_EXP)
    
    V_act_2 = j0_catodo * ((CO2 / CO2_2) * (np.exp(exp_arg_2_pos) - np.exp(exp_arg_2_neg))) - J_catodo
    
    return np.array([V_act_1, V_act_2])