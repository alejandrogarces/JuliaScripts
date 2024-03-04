using LinearAlgebra
using DataFrames
using CSV
using SparseArrays
println("---------------------------------------")
println("Leyendo datos de entrada")
nodos = DataFrame(CSV.File("nodos.csv"))
lineas = DataFrame(CSV.File("lineas.csv"))
num_lineas = nrow(lineas)
num_nodos = nrow(nodos)

# Generar matriz B para flujo DC
B = zeros(num_nodos,num_nodos)
for lin in eachrow(lineas)
    n1 = lin.From
    n2 = lin.To
    y_lin  = 1/lin.x_pu
    B[n1,n1] +=  y_lin
    B[n1,n2] += -y_lin
    B[n2,n1] += -y_lin
    B[n2,n2] +=  y_lin  
end
# eliminar el slack, hacerla dispersa y solucionar el flujo de carga basico
B = sparse(B[2:num_nodos,2:num_nodos])
P = zeros(num_nodos,1)
for nod in eachrow(nodos)
    n = nod.Num
    P[n] += (nod.Generation_MW - nod.Load_MW)/100
end
th = B\P[2:num_nodos]
th = [0;th]  # angulos nodales sin compensacion
# flujos sin compensacion
nodos.Angle = th*180/pi
display(nodos.Angle)
