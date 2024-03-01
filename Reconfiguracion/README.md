## Flujo de carga usando el método de punto fijo

Basado en el siguiente artículo:

*S. Civanlar, J. J. Grainger, H. Yin and S. S. H. Lee, "Distribution feeder reconfiguration for loss reduction," in IEEE Transactions on Power Delivery, vol. 3, no. 3, pp. 1217-1223, July 1988, doi: 10.1109/61.193906.*

Este sistema esta diseñado para reconfiguración de alimentadores primarios.  Falta implementar el algoritmo de reconfiguración.

---
## Características del sistema de prueba

* Potencia base de 100 MW
* Se asume que hay tres nodos slack, correspondientes a las tres subestaciones
* La $Y_\text{bus}$ no tiene en cuenta las lineas cuya columna ''on'' este en cero
* Los datos estan en un dataframe con la siguiente configuracion

|from|to|r_pu|x_pu|p_MW|q_MW|on|
|--|--|--|--|--|--|--|
|1|4|0.075|0.10|2.0|1.6|1|
4|5|0.080|0.11|3.0|0.4|1
4|6|0.090|0.18|2.0|-0.4|1
6|7|0.040|0.04|1.5|1.2|1
2|8|0.110|0.11|4.0|2.7|1
8|9|0.080|0.11|5.0|1.8|1
8|10|0.110|0.11|1.0|0.9|1
9|11|0.110|0.11|0.6|-0.5|1
9|12|0.080|0.11|4.5|-1.7|1
3|13|0.110|0.11|1.0|0.9|1
13|14|0.090|0.12|1.0|-1.1|1
13|15|0.080|0.11|1.0|0.9|1
15|16|0.040|0.04|2.1|-0.8|1
5|11|0.040|0.04|0.0|0.0|0
10|14|0.040|0.04|0.0|0.0|0
7|16|0.090|0.12|0.0|0.0|0


## Funciones
### calcular_ybus
    function calcular_ybus(datos,num_nodos,num_lineas)
Calcula la matriz de admitancia nodal, a partir de los datos de entrada, el número de nodos y el número de lineas.  Regresa una matriz compleja de tamaño $\text{num}_\text{nodos}\times\text{num}_\text{nodos}$
### calcular flujo de carga
    function calcular_flujo_carga(Ybus,S,N,num_nodos,slack)
Método de punto fijo para el flujo de carga $v=T(v)$.  Los datos de entrada son la $Y_\text{bus}$, las potencias nodales, el conjunto de nodos, el numero de nodos y los nodos slack.

### main function
    function main()
Funcion principal.    



