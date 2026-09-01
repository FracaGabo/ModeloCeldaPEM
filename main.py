import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from presion_sat import presion_sat
from Difusividad import Difusividad
from funcion_lambda import funcion_lambda
from resolver_balances import resolver_balances
from sobrevoltaje_activacion import sobrevoltaje_activacion


#CONSTANTES Y PROPIEDADES
R = 8.314 #J/mol K -- Cte Gases Ideales
F = 96485 #C/mol -- Cte de Faraday

PM = (0.002, 0.032, 0.018, 0.028) #Pesos Moleculares -- H2, O2, H2O, N2

parametros = {
    # Parámetros electroquímicos
    "alpha_anodo": 1.0,
    "alpha_catodo": 1.0,
    "j0_anodo": 2700,       # A/m2
    "j0_catodo": 1,         # A/m2
    "J_limite": 23000,      # A/m2

    # Geometría de la celda
    "A_anodo": 2.5e-3,      # m2
    "A_catodo": 2.5e-3,     # m2
    "V_anodo": 0.005,       # m3
    "V_catodo": 0.01,       # m3
    "e_membrana": 50e-4,    # cm
    "rho_membrana": 1980,   # kg/m3
    "PM_membrana": 1100,    # kg/mol
    "n": 110,

    # Condiciones de operación
    "T_operacion": 348.15,  # K
    "exH2": 1.5,
    "exO2": 2.0,
    "P_anodo": 1.3,
    "P_catodo": 2,
    "RH_1": 0.5,
    "RH_2": 0.5,
}

T_operacion = parametros["T_operacion"]
exH2 = parametros["exH2"]
exO2 = parametros["exO2"]

I = 0.8*1e4*parametros["A_anodo"]; #A -- PUEDE VARIAR

#-------------------------------------------------------------------------------------
#CORRIENTE DE HIDROGENO

nH2 = exH2*(I/(2*F))
T_1 = T_operacion
P_1 = parametros["P_anodo"]
Psat_1 = presion_sat(T_1,1)
RH_1 = parametros["RH_1"]

yH2O_1 = (RH_1*Psat_1)/(P_1)
yH2_1 = 1 - yH2O_1

PH2_1 = P_1*yH2_1
PH2O_1 = P_1*yH2O_1

nH2O_1 = nH2/yH2_1 - nH2
CH2_1 = (PH2_1*1e5)/(R*T_1)
CH2O_1 = (PH2O_1*1e5)/(R*T_1)

PM_1 = yH2_1*PM[0]+ yH2O_1*PM[2]


#rho_1 = ((P_1*1e5)*PM_1)/(R*T_1)

q_1 = (nH2 + nH2O_1)*((R*T_1)/(P_1*1e5))


#CORRIENTE DE AIRE ENTRADA
n_O2 = exO2*(I/(4*F))
T_2 = T_operacion
P_2 = parametros["P_catodo"]
RH_2 = parametros["RH_2"]

Psat_2 = presion_sat(T_2,1)

yH2O_2 = (RH_2*Psat_2)/P_2

yair = 1-yH2O_2
yO2_2 = 0.21*yair
yN2_2 = 0.79*yair

PH2O_2 = yH2O_2*P_2
PO2_2 = yO2_2*P_2
PN2_2 = yN2_2*P_2

CH2O_2 = (PH2O_2*1e5)/(R*T_2)
CO2_2 = (PO2_2*1e5)/(R*T_2)
CN2_2 = (PN2_2*1e5)/(R*T_2)

nH2O_2 = (n_O2/yO2_2)*yH2O_2
nair_2 = (n_O2/yO2_2) - nH2O_2
nN2_2 = nair_2 - n_O2

PM_2 = yH2O_2*PM[2] + yO2_2*PM[1] + yN2_2*PM[3]

q_2 = (n_O2 + nH2O_2 + nN2_2)*((R*T_2)/(P_2*1e5))

#rho_2 = ((P_2*1e5)*PM_2)/(R*T_2)

entradas = {
    "q1": q_1,
    "q2": q_2,
    "CH2_1": CH2_1,
    "CH2O_1": CH2O_1,
    
    "CO2_2": CO2_2,
    "CH2O_2": CH2O_2,
    "CN2_2": CN2_2,
    
    "I": I,
}


C_inicial = np.array([
    CH2_1,
    CH2O_1,
    CO2_2,
    CH2O_2,
    CN2_2
])

t, concentraciones, FM = resolver_balances(
    parametros,
    entradas,
    C_inicial,
    0,
    1800,
    dt=1
)


# CÁLCULO DE VOLTAJES
rango = len(t)

# PRESIONES PARCIALES EN BAR
P_parcial = np.zeros((rango, 5))

P_parcial[:, 0] = concentraciones[:, 0] * R * T_operacion * (1e-5)  # Hidrógeno
P_parcial[:, 1] = concentraciones[:, 2] * R * T_operacion * (1e-5)  # Oxígeno
P_parcial[:, 2] = concentraciones[:, 1] * R * T_operacion * (1e-5)  # Agua ánodo
P_parcial[:, 3] = concentraciones[:, 3] * R * T_operacion * (1e-5)  # Agua cátodo
P_parcial[:, 4] = concentraciones[:, 4] * R * T_operacion * (1e-5)  # Nitrógeno

# Inicializar arrays
E_nernst = np.zeros(rango)
sobrevoltaje = np.zeros((rango, 4))
voltaje = np.zeros(rango)
voltaje_stack = np.zeros(rango)
actividad = np.zeros((rango, 2))
lambda_vals = np.zeros((rango, 3))
c_membrana = np.zeros(rango)
Dw = np.zeros(rango)

# POTENCIAL DE REDUCCIÓN ESTÁNDAR
E_standar = 1.229 - 0.9e-3*(T_operacion - 298.15)

# VOLTAJE DE NERNST
for i in range(rango):
    E_nernst[i] = (E_standar - 0.85e-3*(T_operacion - 298.15) + 
                   4.3085e-5 * T_operacion * 
                   (np.log(P_parcial[i, 0]) + 0.5*np.log(P_parcial[i, 1])))


# ============================================================
# SOBREVOLTAJES DE ACTIVACIÓN
# ============================================================

a2 = np.array([
    parametros["j0_anodo"],
    parametros["j0_catodo"],
    CH2_1,
    CO2_2,
    parametros["alpha_anodo"],
    parametros["alpha_catodo"],
    parametros["A_anodo"],
    parametros["A_catodo"]
])


# Primera estimación solamente para t = 0
iterador = np.array([0.05, 0.5])


for i in range(rango):

    t0 = t[i]

    if t0 < 800:
        It = I

    elif t0 < 1200:
        It = 1.6 * parametros["A_catodo"] * 1e4

    else:
        It = 0.8 * parametros["A_catodo"] * 1e4


    # ========================================================
    # VARIABLES DE DECISIÓN
    # ========================================================

    u2 = np.array([
        concentraciones[i, 0],   # H2
        concentraciones[i, 2],   # O2
        T_operacion,
        It
    ])


    V_act, info, ier, mensaje = fsolve(
        sobrevoltaje_activacion,
        iterador,
        args=(u2, a2),
        full_output=True,
        xtol=1e-12,
        maxfev=10000
    )


    # ========================================================
    # COMPROBAR CONVERGENCIA
    # ========================================================

    if ier == 1:

        sobrevoltaje[i, 0] = V_act[0]
        sobrevoltaje[i, 1] = V_act[1]

        # MUY IMPORTANTE:
        # usar esta solución como iterador del siguiente tiempo
        iterador = V_act.copy()

    else:

        print(
            f"fsolve no convergió en t = {t0:.1f} s | "
            f"CH2 = {u2[0]:.6f} | "
            f"O2 = {u2[1]:.6f} | "
            f"I = {It:.3f} A"
        )

        print("Mensaje:", mensaje)

        # Mantener el último valor válido
        if i > 0:

            sobrevoltaje[i, 0] = sobrevoltaje[i-1, 0]
            sobrevoltaje[i, 1] = sobrevoltaje[i-1, 1]

        else:

            sobrevoltaje[i, 0] = iterador[0]
            sobrevoltaje[i, 1] = iterador[1]
    
    
    
    
    
    
    
    # SOBREVOLTAJE OHMICO
    psat = presion_sat(T_operacion, 1)  # bar
    
    actividad[i, 0] = P_parcial[i, 2] / psat #anodo
    actividad[i, 1] = P_parcial[i, 3] / psat #catodo
    lambda_vals[i, :] = funcion_lambda(actividad[i, :])
    Dw[i] = Difusividad(lambda_vals[i, :], T_operacion)
    
    c_membrana[i] = (0.005139 * lambda_vals[i, 2] - 0.00326) * np.exp(1268 * ((1/303) - (1/T_operacion)))
    
    A_cm = parametros["A_anodo"] * 1e4  # cm²
    R_membrana = parametros["e_membrana"] / (c_membrana[i] * A_cm)
    sobrevoltaje[i, 2] = It * R_membrana  # V
    
    J = It / parametros["A_anodo"]
    
    # SOBREVOLTAJE DE CONCENTRACIÓN
    sobrevoltaje[i, 3] = ((R * T_operacion) / (2 * F)) * np.log(parametros["J_limite"] / (parametros["J_limite"] - J))
    
    # VOLTAJE DE LA CELDA
    voltaje[i] = E_nernst[i] - sobrevoltaje[i, 0] - sobrevoltaje[i, 1] - sobrevoltaje[i, 2] - sobrevoltaje[i, 3]
    
    voltaje_stack[i] = voltaje[i] * parametros["n"]

# ...existing code...

# GRAFICO DE CONCENTRACIONES
fig = plt.figure(figsize=(10, 6))
ax = plt.axes()

plt.plot(t, concentraciones[:, 0], label='H₂', linewidth=2)
plt.plot(t, concentraciones[:, 1], label='H₂O (anodo)', linewidth=2)
plt.plot(t, concentraciones[:, 2], label='O₂', linewidth=2)
plt.plot(t, concentraciones[:, 3], label='H₂O (catodo)', linewidth=2)
plt.plot(t, concentraciones[:, 4], label='N₂', linewidth=2)

plt.xlabel('Tiempo (s)', fontsize=12)
plt.ylabel('Concentración (mol/m³)', fontsize=12)
plt.title('Concentración de especies vs Tiempo', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.xlim(0, max(t))
plt.ylim(0, np.max(concentraciones) * 1.1)

plt.tight_layout()

# GRÁFICO 1: VOLTAJE DE CELDA vs VOLTAJE DE NERNST
fig1, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(t, E_nernst, label='Voltaje de Nernst (reversible)', linewidth=2, color='green')
ax1.plot(t, voltaje, label='Voltaje de la celda', linewidth=2, color='blue')
ax1.plot(t, voltaje_stack, label='Voltaje del stack', linewidth=2, color='red')

ax1.set_xlabel('Tiempo (s)', fontsize=12)
ax1.set_ylabel('Voltaje (V)', fontsize=12)
ax1.set_title('Voltaje de Nernst vs Voltaje de la Celda', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max(t))
ax1.set_ylim(0)

plt.tight_layout()

# GRÁFICO 2: SOBREVOLTAJES
fig2, ax2 = plt.subplots(figsize=(10, 6))

ax2.plot(t, sobrevoltaje[:, 0], label='Sobrevoltaje Ánodo (activación)', linewidth=2, color='orange')
ax2.plot(t, sobrevoltaje[:, 1], label='Sobrevoltaje Cátodo (activación)', linewidth=2, color='purple')
ax2.plot(t, sobrevoltaje[:, 2], label='Sobrevoltaje Ohmico', linewidth=2, color='brown')
ax2.plot(t, sobrevoltaje[:, 3], label='Sobrevoltaje Concentración', linewidth=2, color='pink')

ax2.set_xlabel('Tiempo (s)', fontsize=12)
ax2.set_ylabel('Sobrevoltaje (V)', fontsize=12)
ax2.set_title('Sobrevoltajes en la Celda de Combustible', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max(t))
ax2.set_ylim(0)

plt.tight_layout()

datos_exportar = pd.DataFrame({
    'Tiempo (s)': t,
    'H2 (mol/m³)': concentraciones[:, 0],
    'H2O_anodo (mol/m³)': concentraciones[:, 1],
    'O2 (mol/m³)': concentraciones[:, 2],
    'H2O_catodo (mol/m³)': concentraciones[:, 3],
    'N2 (mol/m³)': concentraciones[:, 4],
    'P_parcial_H2 (bar)': P_parcial[:, 0],
    'P_parcial_O2 (bar)': P_parcial[:, 1],
    'P_parcial_H2O_anodo (bar)': P_parcial[:, 2],
    'P_parcial_H2O_catodo (bar)': P_parcial[:, 3],
    'P_parcial_N2 (bar)': P_parcial[:, 4],
    'Actividad_anodo': actividad[:, 0],
    'Actividad_catodo': actividad[:, 1],
    'Lambda_1': lambda_vals[:, 0],
    'Lambda_2': lambda_vals[:, 1],
    'Lambda_3': lambda_vals[:, 2],
    'Difusividad (m2/s)': Dw,
    'Conductividad_membrana (S/cm)': c_membrana,
    'Voltaje_Nernst (V)': E_nernst,
    'Sobrevoltaje_Anodo (V)': sobrevoltaje[:, 0],
    'Sobrevoltaje_Catodo (V)': sobrevoltaje[:, 1],
    'Sobrevoltaje_Ohmico (V)': sobrevoltaje[:, 2],
    'Sobrevoltaje_Concentracion (V)': sobrevoltaje[:, 3],
    'Voltaje_Celda (V)': voltaje,
    'Voltaje_Stack (V)': voltaje_stack,
})


# Guardar en Excel
ruta_excel = r'd:\Gab\Escritorio\Proyectos_Python\resultados_3.xlsx'
datos_exportar.to_excel(ruta_excel, index=False, sheet_name='Datos')

print(f"✓ Datos exportados exitosamente a: {ruta_excel}")
print(f"✓ Filas: {len(datos_exportar)}")
print(f"✓ Columnas: {len(datos_exportar.columns)}")
print("✓ Simulación completada")


# Mostrar todos los gráficos juntos al final
plt.show()