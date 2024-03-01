# Analisis de estabilidad transitoria usando el modelo clásico

## Resumen

Estabilidad transitoria usando el modelo clásico.  Basado en el sistema de preba del libro

*Paul M. Anderson; A. A. Fouad, "Frontmatter," in Power System Control and Stability , IEEE, 2003, pp.i-xiv, doi: 10.1109/9780470545577.fmatter.*
## Funciones

La mayor parte del codigo esta en el script principal.  Solo se usaron las siguientes funciones:

    function kron_reduction(Y::Matrix{ComplexF64},n::Int64)

Calcula la reducción de Kron

    function Base.show(io::IO, vec_d::Vector{NamedTuple})

imprime el diccionario correspondiente

---
## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
