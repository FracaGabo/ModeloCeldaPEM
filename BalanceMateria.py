import numpy as np 
from scipy.integrate import solve_ivp

#BALANCE MATERIA HIDRÓGEN

def balances_materia(t,Ci, parametros):
    
    CH2, CO2, CH2O_AN, CH2O_CA, CN2 = Ci
    
    R = 8.314
    
    q1 = parametros["q1"]
    q2 = parametros["q2"]
    I = parametros["I"]
    F = parametros["F"]
    FM = parametros["FM"]

    CH2_1 = parametros["CH2_1"]
    CH2O_1 = parametros["CH2O_1"]

    CO2_2 = parametros["CO2_2"]
    CH2O_2 = parametros["CH2O_2"]
    CN2_2 = parametros["CN2_2"]

    V_anodo = parametros["V_anodo"]
    V_catodo = parametros["V_catodo"]
    
    T_operacion = parametros["T_operacion"]
    
    
    #ITERADORES INICIALES
    FM = 1e-4
    iterador = 1
    tolerancia = 1e-8
    max_iteraciones = 1e5
    error = 1e10
    
    while iterador <= max_iteraciones and error > tolerancia:
        
    
        q3 = ((q1*CH2_1-I/(2*F))+(q1*CH2O_1+FM))/(CH2+CH2O_AN)
    
        q4 = (q2*CN2_2 + q2*CO2_2 - (I/(4*F)) + q2*CH2O_2 + (I/(2*F)) - FM)/(CN2 + CO2 +CH2O_CA)
    
        #BALANCE DE MATERIA ANODO
        
        dCH2dt = (q1*CH2_1 - q3*CH2 - I/(2*F))/V_anodo
        dCH2O_ANdt = (q1*CH2O_1 - q3*CH2O_AN + FM)/V_anodo    
        
        #BALANCE DE MATERIA CATODO
        
        dCO2dt = (q2*CO2_2 - q4*CO2 - I/(4*F))/V_catodo
        dCH2O_CAdt = (q2*CH2O_2 -q4*CH2O_CA + I/(2*F) - FM)/V_catodo
        dCN2dt = (q2*CN2_2 - q4*CN2)/V_catodo
    
        presion_saturacion = presion_sat(T_operacion,1)*1e5 #Pa
        
        actividad_agua = [((R*T_operacion*CH2O_AN)/presion_saturacion), ((R*T_operacion*CH2O_CA)/presion_saturacion)]
        
        funcion_lambda = contenido_agua(actividad_agua)
        
        Dw = difusividad_agua(funcion_lambda,T_operacion)
        
        nd = (2.5/22)
        
        CH2O_AN = (rho_membrana/PM_membrana)*funcion_lambda[0]
        CH2O_CA = (rho_membrana/PM_membrana)*funcion_lambda[1]

        FM_calculado = nd*(I/F) - ((A_catodo*Dw*(CH2O_AN-CH2O_CA)))/(e_membrana)
        
        error = abs(FM-FM_calculado)
        
        FM = FM_calculado
        iterador += 1
        
    
    
    
    return(dCH2dt, dCO2dt, dCH2O_ANdt, dCH2O_CAdt, dCN2dt)