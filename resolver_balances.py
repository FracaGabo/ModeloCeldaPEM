import numpy as np
from scipy.integrate import solve_ivp


def resolver_balances(parametros, entradas, C_inicial, ti, tf, dt=1):

    # =========================================================
    # PARÁMETROS
    # =========================================================

    F = parametros["F"]
    R = parametros["R"]

    V_anodo = parametros["V_anodo"]
    V_catodo = parametros["V_catodo"]

    A_catodo = parametros["A_catodo"]
    e_membrana = parametros["e_membrana"]

    rho_membrana = parametros["rho_membrana"]
    PM_membrana = parametros["PM_membrana"]


    # =========================================================
    # ENTRADAS
    # =========================================================

    q1 = entradas["q1"]
    q2 = entradas["q2"]

    CH2_1 = entradas["CH2_1"]
    CH2O_1 = entradas["CH2O_1"]

    CO2_2 = entradas["CO2_2"]
    CH2O_2 = entradas["CH2O_2"]
    CN2_2 = entradas["CN2_2"]

    T_operacion = entradas["T_operacion"]


    # =========================================================
    # VECTOR DE TIEMPO
    # =========================================================

    tiempos = np.arange(ti, tf + dt, dt)


    # =========================================================
    # MATRICES PARA GUARDAR RESULTADOS
    # =========================================================

    concentraciones = np.zeros((len(tiempos), 5))

    FM_resultados = np.zeros(len(tiempos))


    # =========================================================
    # CONDICIONES INICIALES
    # =========================================================

    concentraciones[0, :] = C_inicial


    # =========================================================
    # VALOR INICIAL DEL ITERADOR FM
    # =========================================================

    FM_anterior = 1e-4

    FM_resultados[0] = FM_anterior


    # =========================================================
    # PARÁMETROS DE CONVERGENCIA
    # =========================================================

    tolerancia = 1e-5

    max_iteraciones = 1000


    # =========================================================
    # CICLO TEMPORAL
    # =========================================================

    for k in range(len(tiempos) - 1):

        # Tiempo inicial y final del intervalo
        t0 = tiempos[k]
        t1 = tiempos[k + 1]

        # Concentraciones al inicio del intervalo
        C0 = concentraciones[k, :].copy()

        # Primera estimación de FM
        FM_iter = FM_anterior

        convergio = False


        # =====================================================
        # ITERACIÓN PARA HACER CONVERGER FM
        # =====================================================

        for iteracion in range(max_iteraciones):


            # =================================================
            # CORRIENTE EN ESTE INTERVALO
            # =================================================

            if t0 < 800:

                I = entradas["I"]

            elif t0 < 1200:

                I = 1.6 * A_catodo * 1e4

            else:

                I = 0.8 * A_catodo * 1e4


            # =================================================
            # FUNCIÓN DE BALANCES DE MATERIA
            # =================================================

            def balances(t, C):

                CH2, CH2O_AN, CO2, CH2O_CA, CN2 = C


                # CAUDALES DE SALIDA

                q3 = (q1*CH2_1 - I/(2*F) + q1*CH2O_1 - FM_iter) / (CH2O_AN + CH2)

                q4 = (q2*CN2_2 + q2*CO2_2 - I/(4*F) + q2*CH2O_2 + I/(2*F) + FM_iter) / (CN2 + CO2 + CH2O_CA)


                # =============================================
                # BALANCE DE MATERIA ÁNODO
                # =============================================

                # HIDRÓGENO

                dCH2dt = (q1/V_anodo)*CH2_1 - (q3/V_anodo)*CH2 - I/(2*F*V_anodo)


                # AGUA ÁNODO

                dCH2O_ANdt = (q1/V_anodo)*CH2O_1 - (q3/V_anodo)*CH2O_AN - FM_iter/V_anodo


                # =============================================
                # BALANCE DE MATERIA CÁTODO
                # =============================================

                # OXÍGENO

                dCO2dt = (q2/V_catodo)*CO2_2 - (q4/V_catodo)*CO2 - I/(4*F*V_catodo)


                # AGUA CÁTODO

                dCH2O_CAdt = (q2/V_catodo)*CH2O_2 - (q4/V_catodo)*CH2O_CA + I/(2*F*V_catodo) + FM_iter/V_catodo


                # NITRÓGENO

                dCN2dt = (q2/V_catodo)*CN2_2 - (q4/V_catodo)*CN2


                return [
                    dCH2dt,
                    dCH2O_ANdt,
                    dCO2dt,
                    dCH2O_CAdt,
                    dCN2dt
                ]


            # =================================================
            # RESOLUCIÓN DE LAS EDO PARA ESTE INTERVALO
            # =================================================

            solucion = solve_ivp(
                balances,
                (t0, t1),
                C0,
                method="BDF",
                rtol=1e-6,
                atol=1e-8
            )


            # Verificar que solve_ivp pudo resolver el intervalo

            if not solucion.success:

                raise RuntimeError(
                    f"Error en solve_ivp en t = {t0} s: "
                    f"{solucion.message}"
                )


            # =================================================
            # CONCENTRACIONES AL FINAL DEL INTERVALO
            # =================================================

            C_nueva = solucion.y[:, -1]


            CH2 = C_nueva[0]

            CH2O_AN = C_nueva[1]

            CO2 = C_nueva[2]

            CH2O_CA = C_nueva[3]

            CN2 = C_nueva[4]


            # =================================================
            # PRESIÓN DE SATURACIÓN
            # =================================================

            psat_Pa = presion_sat(T_operacion, 1) * 1e5


            # =================================================
            # ACTIVIDAD DEL AGUA
            # =================================================

            actividad = [
                R*T_operacion*CH2O_AN/psat_Pa,
                R*T_operacion*CH2O_CA/psat_Pa
            ]


            # =================================================
            # CONTENIDO DE AGUA DE LA MEMBRANA
            # =================================================

            lambda_agua = contenido_agua(actividad)


            # =================================================
            # DIFUSIVIDAD DEL AGUA
            # =================================================

            Dw = difusividad(lambda_agua, T_operacion)


            # =================================================
            # CONCENTRACIÓN DE AGUA SEGÚN LA MEMBRANA
            # =================================================

            CH2O_AN_mem = (rho_membrana/PM_membrana) * lambda_agua[0]

            CH2O_CA_mem = (rho_membrana/PM_membrana) * lambda_agua[1]


            # =================================================
            # ARRASTRE ELECTROOSMÓTICO
            # =================================================

            nd = 2.5/22


            # =================================================
            # NUEVO FM
            # =================================================

            FM_nuevo = nd*(I/F) - A_catodo*Dw*(CH2O_CA_mem - CH2O_AN_mem)/e_membrana


            # =================================================
            # ERROR ENTRE FM ANTERIOR Y FM NUEVO
            # =================================================

            error = abs(FM_nuevo - FM_iter)


            # =================================================
            # COMPROBAR CONVERGENCIA
            # =================================================

            if error < tolerancia:

                FM_iter = FM_nuevo

                convergio = True

                break


            # Si no converge, FM_nuevo pasa a ser
            # el FM de la siguiente iteración

            FM_iter = FM_nuevo


        # =====================================================
        # VERIFICAR CONVERGENCIA
        # =====================================================

        if not convergio:

            print(
                f"Advertencia: FM no convergió "
                f"en t = {t1} s"
            )


        # =====================================================
        # GUARDAR CONCENTRACIONES CONVERGIDAS
        # =====================================================

        concentraciones[k + 1, :] = C_nueva


        # =====================================================
        # GUARDAR FM CONVERGIDO
        # =====================================================

        FM_resultados[k + 1] = FM_iter


        # =====================================================
        # USAR FM CONVERGIDO COMO ESTIMACIÓN DEL SIGUIENTE PASO
        # =====================================================

        FM_anterior = FM_iter


    # =========================================================
    # RESULTADOS
    # =========================================================

    return tiempos, concentraciones, FM_resultados