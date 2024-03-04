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
    B[n2,n2] += y_lin  
end
# eliminar el slack, hacerla dispersa y solucionar el flujo de carga basico
B = sparse(B[2:num_nodos,2:num_nodos])
P = zeros(num_nodos,1)
for nod in eachrow(nodos)
    n = nod.Num
    P[n] += (nod.Generation_MW - nod.Load_MW)/100
end
th = B\P[2:num_nodos]
th = [0;th]  # angulos nodales

# flujos y perdidas aproximadas
p_f = zeros(num_lineas)
p_loss = zeros(num_lineas)
for (k,lin) in enumerate(eachrow(lineas))
    th_km = th[lin.From]-th[lin.To]
    g_km = lin.r_pu/(lin.r_pu^2+lin.x_pu^2)
    p_f[k] = th_km/lin.x_pu
    p_loss[k] = g_km*(th_km^2)
end

# resultados
p_loss_total = sum(p_loss)*100
println("Perdidas aprox:",p_loss_total)
println("Potencia slack:",p_loss_total-sum(P[2:num_nodos])*100)
nodos.Angle = th*180/pi
display(nodos.Angle)
