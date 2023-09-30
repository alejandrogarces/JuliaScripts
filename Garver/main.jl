using LinearAlgebra
using DataFrames
using CSV

lines  = DataFrame(CSV.File("data_lines.csv"))
nodes  = DataFrame(CSV.File("data_nodes.csv"))
num_nodes = maximum(nodes.node)
num_lines = nrow(lines)

using JuMP
using HiGHS
model = Model(HiGHS.Optimizer)
@variable(model, -pi/2 <= theta[1:num_nodes] <= pi/2 )
@variable(model, p_lin[1:num_lines])
@variable(model,z[1:num_lines],Bin)
A = zeros(num_nodes,num_lines)
for k in 1:num_lines
    i = lines.from[k]
    j = lines.to[k]
    @constraint(model, p_lin[k] <=  z[k]*lines.pmax_MW[k]/100)
    @constraint(model, p_lin[k] >= -z[k]*lines.pmax_MW[k]/100)
    @constraint(model, p_lin[k]*lines.x_pu[k] - (theta[i]-theta[j]) <=  (1-z[k])*pi)
    @constraint(model, p_lin[k]*lines.x_pu[k] - (theta[i]-theta[j]) >= -(1-z[k])*pi)
    A[i,k] =  1
    A[j,k] = -1
end
@constraint(model, theta[1] == 0.0)
Pn = (nodes.gen_MW[:]-nodes.load_MW[:])/100
@constraint(model,Pn .== A*p_lin)
@objective(model,Min,sum(lines.cost.*z))
optimize!(model)
println("Costo total = ", objective_value(model))
println("Lines to connect ")
for k = 1:num_lines
    if value(z[k])>0
        i = lines.from[k]
        j = lines.to[k]
        if lines.cost[k] > 0
          printstyled("line ",k,": ",i,"-->",j, "  ", round(value(p_lin[k])*100,digits=2)," MW\n", color=:blue)
        else
            println("line ",k,": ",i,"-->",j, "  ", round(value(p_lin[k])*100,digits=2)," MW")
        end
    end
end
