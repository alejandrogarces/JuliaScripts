# Estimación de estados usando un método de punto fijo para redes DC

La función principal toma una microrred DC y calcula los voltajes usando un flujo de carga convencional.  Posteriormente, introduce ruido en las mediciones de voltaje y potencia nodal.  El algoritmo de optimización de punto fijo es usado para estimar el estado del sistema a partir de estas medidas distorcionadas. Se usa un modelo de optimizacion en JuMP para comparar.  Igualmente, se desarrolla un algoritmo de estimacion de la constante Lipschitz con el objetivo de justificar de que en efecto T es un mapa de contraccion.

## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
