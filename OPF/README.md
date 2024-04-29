# Flujo de carga óptimo usando linealización de Wirtinger

Resuelve el problema de flujo de carga óptimo usando una linealización sobre los numero complejos.  El modelo de optimización es el siguiente

$\min V_r^\top G V_r + V_i^\top G V_i$

sujeto a:

$v_k^* = y_{km} (v_k^*+v_m-1) $ 

$\left\|v_{k}-1\right\|_1 \leq 0.1$

$0 \leq p_{g} \leq p_{max}$ 

$0 \leq q_{g} \leq q_{max}$ 


## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
