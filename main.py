#CONSTANTES Y PROPIEDADES
R = 8.314 #J/mol K -- Cte Gases Ideales
F = 96485 #C/mol -- Cte de Faraday

PM = (0.002, 0.032, 0.018, 0.028) #Pesos Moleculares -- H2, O2, H2O, N2

alpha_anodo = 1; # adimensional
alpha_catodo = 1; # adiomensional

j0_anodo = 2700; #A/m2
j0_catodo = 1; #A/m2

J_limite = 23000; #A/m2
 
#GEOMETRÍA DE LA CELDA
A_anodo = 2.5*1e-3 #m2
A_catodo = 2.5*1e-3 #m2
V_anodo = 0.005 #m3
V_catodo = 0.01 #m3
e_membrana = 50*1e-4 #cm
rho_membrana = 1980 #kg/m3
PM_membrana = 1100 #kg/mol
n = 110;

#CONDICIONES DE OPERACIÓN

T_operacion = 348.15; #K
exH2 = 1.5;
exO2 = 2;

I = 5; #A -- PUEDE VARIAR

