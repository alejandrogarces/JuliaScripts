using LinearAlgebra
using DataFrames
# conductance matrix
function Conductance(Lines,num_lines,num_nodes)
    G = zeros(num_nodes,num_nodes)
    for line in eachrow(Lines)
        gkm = 1/line.rkm
        G[line.k,line.k] +=  gkm
        G[line.k,line.m] += -gkm
        G[line.m,line.k] += -gkm
        G[line.m,line.m] +=  gkm
    end
    return G
end

# main function
# input data
Lines = DataFrame( k = [1,2,3],
                   m = [2,3,4],
                   rkm = [0.125,0.5,0.250])
num_lines = size(Lines,1)
num_nodes = maximum([Lines.k;Lines.m])
# per unit representation
v_base = 100; p_base = 100; r_base = v_base^2/p_base;    
Lines.rkm = Lines.rkm/r_base;
G = Conductance(Lines,num_lines,num_nodes)

display(G)
N = 1:num_nodes
Nodes = DataFrame(p = [0.0,200.0,150.0,250.0])
Nodes.p = Nodes.p/p_base

using JuMP
using Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, v[1:num_nodes] >= 0)
@constraint(model,v[1]==1)
p_n = Nodes.p  # todos son generadores
f_loss = 1/2*v'*G*v
@NLobjective(model, Min, f_loss-sum(p_n[k]*log(v[k]) for k in N))
optimize!(model)
v_n = [value(v[k]) for k in N]
p_n[1] = v_n[1]*G[1,:]'*v_n
display(v_n)
H = G+Diagonal(p_n./v_n.^2)
E = eigen(H[2:num_nodes,2:num_nodes])
if minimum(E.values) >= 0
   print("La Hesiana es semidefinida positiva")
end
