# Planeación de la transmision

Modelo de programación lineal para la planeación de la transmisión. 
Basado en el sistema Garver.  Los datos de entrada son información
de las líneas incluyendo costos.  Las líneas existentes tienen costo cero
El modelo es de programación lineal entera mixta, a saber:

Variables de decision:
* $p_k$ flujo por la linea $k=(i,j)$
* $z_k$ variable binaria indicando si la línea existe
* $\theta_i$ angulo en el nodo $i$

## Modelo de optimización:

$\min \sum_k c_kz_k$

sujeto a

$|\theta_{ij}|\leq \pi/2$

$|p_k|\leq z_kp_k^\text{max}$

$|x_k p_k - \theta_{ij}| \leq (1-z_k)\pi$

$[A_{ik}][p_k] = [s_i]$

La última restriccion es matricial, en donde $A$ es la matriz de incidencia.  Note que cuando $z_k=1$ se activa tanto la restriccion de potencia maxima, como la ecuación de flujo DC.  En caso contrario, $p_k=0$ y $\theta_{km}$ queda libre, aunque desde luego no puede tomar un valor por fuera del rango $[-\pi/2,\pi/2]$.

## Datos de entrada

Los datos de entrada son dos archivos data_lines y data_nodes
El primero contiene los datos de las líneas.  Como pueden haber varias lineas en un mismo corredor, se repite el dato las veces que se permitan las líneas en dicho corredor.  Los datos nodales son basicamente potencia nodal.  Se divide por una base de 100 MW



