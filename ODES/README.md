# Comparación entre diferentes métodos de solucion de ecuaciones diferenciales

$f(x) = (J(x)-R(x))\nabla H(x) + G(x) u(x)$

## Método de Euler explícito

$x^{+} = x + \tau f(x)$

## Método de Euler implicito

$x^{+} = x + \tau f(x^{+})$

## Método de Runge-Kuta explicito

$f_1 = f(x)$

$f_2 = f(x+\tau f_1/2)$

$f_3 = f(x+\tau f_2/2)$

$f_4 = f(x+\tau f_3)$

$x^{+} = x + \frac{\tau}{6} (f_1 + 2f_2 + 2f_3 + f_4)$


## Integración de Verlet
FALTA

## Método del salto de rana
FALTA


---
## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
