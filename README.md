# Modelo dinamico PEMFC

Simula balances de materia, transporte de agua en la membrana y voltaje de una
celda de combustible PEM sometida a cambios de corriente.

## Unidades

El nucleo del modelo usa SI. Las presiones de operacion se ingresan en bar y la
conductividad de membrana se conserva en S/cm porque asi esta definida su
correlacion. El espesor se almacena una sola vez en metros y se convierte solo
al calcular la resistencia electrica.

## Ejecucion

```powershell
python main.py
```

El programa exporta `resultados_3.xlsx` antes de mostrar los graficos.

## Pruebas

```powershell
python -m unittest discover -s tests -v
```

Los parametros se encuentran en `config.py`. La opcion
`alimentacion_sigue_corriente` permite elegir entre caudales fijos y caudales
que mantienen los excesos de reactivos durante los escalones de corriente.
