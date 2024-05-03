# Comparación entre diferentes métodos de solucion de ecuaciones diferenciales

## Tipo de dato 

Estructura de sistema Hamiltoniano controlado por puertos

$\dot{x} = (J(x)-R(x))\nabla H(x) + G(x) u(x)$
$y = G(x)^\top \nabla H(x)$

    struct pHsystem
        name::String
        J::Function
        R::Function
        dH::Function
        H::Function
        G::Function
        u::Function
        nx :: Int
        nu :: Int 
    end

## Paralel pH systems

Conecta dos sistemas pH en paralelo.  El resultado es un nuevo sistema pH

FALTA:

## Feedback Ph systems

Conecta dos sistema pH realimentados.  El resultado es un nuevo sistema pH


## Pre y post-multiplicacion de sistemas pH

## Método de Euler explícito

Encuentra una trayectoria usando el método de Euler Explicito:

$x^{+} = x + \tau f(x)$

en donde:

$f(x) = (J(x)-R(x))\nabla H(x) + G(x) u(x)$

La estructura de la función es la siguiente:

    function ExplicitEuler(S,xini,dt,nt)

En donde S es un sistema hamiltoniano controlado por puertos, xini es un punto inicial de x, dt es el paso y nt el número de pasos. La salidas son las siguientes:

    return tode,xode,uode,yode,Hode



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

