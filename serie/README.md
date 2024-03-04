# Compensación serie

## Resumen

Minimiza la cargabilidad en unas lineas asociada a un área de operación, considerando compensación serie.  Se usa el sistema de prueba de 30 nodos

## Modelo

Funcion objetivo:  minimizar la cargabilidad de las líneas asociadas al área de operación A:

$\min \sum_{km \in A} (p_{km}/s_{km})^2$

Restricciones:  las convencionales asociadas al flujo de carga:

$-\theta_{km}^{\text{max}} \leq \theta_{km} \leq \theta_{km}^{\text{max}}$

$\theta_{km} - x_{km}p_{km} = 0$

además, una restriccion asociada a las líneas con compensacion serie:

$|\theta_{km}-x_{km}p_{km}| \leq p_{km}x_{\text{max}}$

---
## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

