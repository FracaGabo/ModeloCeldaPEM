import unittest

import numpy as np

from alimentacion import calcular_entradas
from config import corriente_programada, crear_parametros
from electroquimica import (
    potencial_nernst,
    sobrevoltaje_butler_volmer,
    sobrevoltaje_concentracion,
)
from funcion_lambda import funcion_lambda
from main import ejecutar_simulacion
from presion_sat import presion_sat


class PruebasModelo(unittest.TestCase):
    def setUp(self):
        self.parametros = crear_parametros()

    def test_presion_saturacion_agua_75_c(self):
        self.assertAlmostEqual(presion_sat(348.15), 0.385, delta=0.015)

    def test_lambda_en_actividad_cero_es_finito(self):
        valores = funcion_lambda([0.0, 0.0])
        np.testing.assert_allclose(valores, [0.043, 0.043, 0.043])

    def test_perfil_corriente(self):
        base = 20.0
        self.assertEqual(corriente_programada(799, self.parametros, base), 20.0)
        self.assertEqual(corriente_programada(800, self.parametros, base), 40.0)
        self.assertEqual(corriente_programada(1200, self.parametros, base), 20.0)

    def test_nernst_es_finito(self):
        valor = potencial_nernst(348.15, 1.0, 0.21)
        self.assertTrue(np.isfinite(valor))
        self.assertGreater(valor, 1.0)

    def test_alpha_igual_a_uno_es_valido(self):
        eta = sobrevoltaje_butler_volmer(
            corriente=20.0,
            area=2.5e-3,
            j0=2700.0,
            concentracion=10.0,
            concentracion_ref=10.0,
            alpha=1.0,
            T=348.15,
        )
        self.assertGreater(eta, 0.0)

    def test_sobrevoltaje_concentracion_formula_explicita(self):
        densidad = 8000.0
        limite = 23000.0
        calculado = sobrevoltaje_concentracion(348.15, densidad, limite)
        esperado = (8.314462618 * 348.15 / (2 * 96485.33212)) * np.log(
            limite / (limite - densidad)
        )
        self.assertAlmostEqual(calculado, esperado, places=12)

    def test_entradas_positivas(self):
        entradas = calcular_entradas(self.parametros, 20.0)
        for nombre in ("q1", "q2", "CH2_1", "CH2O_1", "CO2_2", "CH2O_2", "CN2_2"):
            self.assertGreater(entradas[nombre], 0)

    def test_simulacion_corta(self):
        resultados = ejecutar_simulacion(self.parametros, tf=2, dt=1)
        self.assertEqual(resultados["concentraciones"].shape, (3, 5))
        self.assertTrue(np.all(resultados["concentraciones"] > 0))
        self.assertTrue(np.all(np.isfinite(resultados["voltaje"])))

    def test_presiones_totales_se_mantienen(self):
        resultados = ejecutar_simulacion(self.parametros, tf=2, dt=1)
        presion_anodo = resultados["P_parcial"][:, 0] + resultados["P_parcial"][:, 1]
        presion_catodo = np.sum(resultados["P_parcial"][:, 2:5], axis=1)
        np.testing.assert_allclose(presion_anodo, self.parametros["P_anodo"], rtol=1e-8)
        np.testing.assert_allclose(presion_catodo, self.parametros["P_catodo"], rtol=1e-8)


if __name__ == "__main__":
    unittest.main()
