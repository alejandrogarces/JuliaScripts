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
lin_compensacion = [3,18,28]
x_comp_max = 0.5     # maxima compensacion capacitiva
lineas_area =      [  1,   2,   3,   4,   5,   6,   7,   8,   9, 10]
lineas_capacidad = [1.2, 0.8, 0.7, 0.7, 0.8, 0.7, 1.0, 0.5, 0.5, 0.5]

# renovables % aproximadamente el 10% de la demanda
nodos.Generation_MW[10:num_nodos] = nodos.Load_MW[10:num_nodos].*rand(num_nodos-9); 

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
th = [0;th]  # angulos nodales sin compensacion
# flujos sin compensacion
p_f = zeros(num_lineas)
for k =  1:num_lineas
    p_f[k] = (th[lineas.From[k]]-th[lineas.To[k]])/lineas.x_pu[k]
end
obj_sin_compensacion = sum((p_f[lineas_area]./lineas_capacidad[lineas_area]).^2)


# Optimizacion de la compensacion serie
theta_max = pi/6
using JuMP
using Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, theta[1:num_nodos])
@variable(model, p_flujo[1:num_lineas])
@variable(model, z_flujo[1:num_lineas])
eq_nodo = Array{GenericAffExpr{Float64, VariableRef}}(undef,1,num_nodos)
for k = 1:num_nodos
    eq_nodo[k] = AffExpr(0.0)
    @constraint(model,theta[k] <= theta_max)
    @constraint(model,theta[k] >= -theta_max)
end

for k = 1:num_lineas
    n1 = lineas.From[k]
    n2 = lineas.To[k]
    x  = lineas.x_pu[k]
    @constraint(model,z_flujo[k] >=  0)
    @constraint(model,p_flujo[k] <=  z_flujo[k])
    @constraint(model,p_flujo[k] >= -z_flujo[k])
    if k in lin_compensacion
       println("línea ",k," : ",abs(p_f[k]))
       @constraint(model, theta[n1]-theta[n2]-x*p_flujo[k] <= z_flujo[k]*x_comp_max)
       @constraint(model, theta[n1]-theta[n2]-x*p_flujo[k] >= -z_flujo[k]*x_comp_max)
    else
       @constraint(model,x*p_flujo[k] == theta[n1]-theta[n2])
    end
    add_to_expression!(eq_nodo[n1], 1,p_flujo[k])
    add_to_expression!(eq_nodo[n2],-1,p_flujo[k])  
end
@variable(model,p_slack)
@constraint(model,p_slack == eq_nodo[1])
for k = 2:num_nodos
    @constraint(model,P[k] == eq_nodo[k])
end

@objective(model,Min,sum((z_flujo[k]/lineas_capacidad[k])^2 for k in lineas_area))
set_silent(model)
optimize!(model)
println(termination_status(model))

println("Flujos por las lineas")
p_c = zeros(num_lineas)
for k in 1:num_lineas
    p_c[k] = value(p_flujo[k])
end
res = DataFrame()
res.flujos_sin_compensacion = p_f
res.flujos_con_compensacion = p_c
display(res)
println("compensaciones")
for k in lin_compensacion
    th_km = value(theta[lineas.From[k]])-value(theta[lineas.To[k]])
    println("linea ",k,",\t x_km+x_c:  ", round(lineas.x_pu[k]-th_km/p_c[k],digits=4))
    lineas.x_pu[k] = th_km/p_c[k] 
end
println("Cargabilidad sin compensacion, :", obj_sin_compensacion)
println("Cargabilidad con compensacion, :", objective_value(model))

using Plots
plt = bar(lineas_area,
          res.flujos_sin_compensacion[lineas_area], 
          color=:gray,
          opacity=0.6, 
          bar_width=0.45) 
plt = bar!(lineas_area.+0.5,
           res.flujos_con_compensacion[lineas_area], 
           color=:royalblue1,
           opacity=0.6,
           bar_width=0.45) 
display(plt)
