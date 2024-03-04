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
lin_compensacion = [5,7]
compensacion = [(5,0.5),(7,0.5)]     # líneas con compensacion y maxima valor capacitivo
lineas_area =      [  1,   2,   3,   4,   5,   6,   7,   8,   9, 10]
lineas_capacidad = [2.0, 1.0, 1.0, 1.0, 1.5, 1.0, 1.0, 1.0, 0.5, 0.5]
# contingencias
#lineas.x_pu[2] = 1E6;

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
for k in lineas_area
    p_f[k] = (th[lineas.From[k]]-th[lineas.To[k]])/lineas.x_pu[k]
end


# Optimizacion de la compensacion serie
theta_max = pi/6
using JuMP
using HiGHS
model = Model(HiGHS.Optimizer)
@variable(model, theta[1:num_nodos])
@variable(model, p_flujo[1:num_lineas])
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
    if k in lin_compensacion
        println("linea ",k)
#        @constraint(model, theta[n1]-theta[n2]-p_flujo[k] >= 0)
#        @constraint(model, theta[n1]-theta[n2]-p_flujo[k] <= -p_flujo[k]*0.05) 
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

@objective(model,Min,sum((p_flujo[lineas_area]./lineas_capacidad[lineas_area]).^2))
set_silent(model)
optimize!(model)
println(termination_status(model))

# imprimir resultados
using Plots
pt1 = bar(1:10,abs.(p_f))
pt1 = bar!(1:10,abs.(p_f))

display(pt1)
println("Flujos por las lineas")
println("_____________________________")
res = DataFrame()
res.From = lineas.From[lineas_area]
res.To = lineas.To[lineas_area]
n_a = length(lineas_area)
p_f = zeros(n_a)
p_c = zeros(n_a)
for k in lineas_area
    p_c[k] = value(p_flujo[k])
end
res.cargabilidad_sin = round.(100*p_f./lineas_capacidad)
res.cargabilidad_con = round.(100*p_c./lineas_capacidad)
#println(res)
println("suma de potencia nodal:", -sum(P[2:num_nodos]))
println("potencia en el slack:", value(p_slack))