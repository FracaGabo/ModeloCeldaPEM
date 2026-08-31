def presion_sat(T,n):
    
    T_C = T - 273.15
    
    if n == 0:
        A = 8.07131
        B = 1730.63
        C = 233.426
        
        Psat = 10**(A-B/(C+T_C))
        Psat = Psat*133.322
        Psat = Psat*1e-5
        
    else:
        A = -2.1794
        B = 0.02953
        C = -9.1837**(1e-5)
        D = 1.4454**(1e-7)
        
        Psat = 10^(A + B*(T_C) + C*(T_C**2) + D*(T_C**3))
        
        
    return Psat