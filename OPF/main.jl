using CSV
using DataFrames
using LinearAlgebra
using JuMP
using Ipopt

function CalculateYbus(Lines)
    n = nrow(Lines)+1
    Ybus = zeros(n,n)*1im    
    for line in eachrow(Lines)
        nt1 = line.From
        nt2 = line.To
        yL  = 1/(line.x*1im + line.r)
        Ybus[nt1,nt1] = Ybus[nt1,nt1] + yL
        Ybus[nt1,nt2] = Ybus[nt1,nt2] - yL
        Ybus[nt2,nt1] = Ybus[nt2,nt1] - yL
        Ybus[nt2,nt2] = Ybus[nt2,nt2] + yL   
    end
    return Ybus
end


function LinearOPF(Lines,Nodes)
    n = nrow(Lines)+1
    Ybus = CalculateYbus(Lines)
    Gbus = real(Ybus)
    model = Model(Ipopt.Optimizer)
    @variable(model,Pg[1:n])
    @variable(model,Qg[1:n])
    @variable(model,Vreal[1:n])
    @variable(model,Vimag[1:n])
    Vn = Vreal + Vimag*1im
    Sn = Pg+Qg*1im-Nodes.Pload - Nodes.Qload*1im
    for k = 1:n
        @constraint(model,conj(Sn[k]) == sum(Ybus[k,m]*(conj(Vn[k])+Vn[m]-1) for m in 1:n))
        @constraint(model,0<= Pg[k] <= Nodes.Pmax[k])
        @constraint(model,-Nodes.Qmax[k]<= Qg[k] <= Nodes.Qmax[k])
        @constraint(model, 0.9 <= Vreal[k] <= 1.1)
        @constraint(model,-0.1 <= Vimag[k] <= 0.1)
    end
    @constraint(model,Vn[1]==1)
    @objective(model,Min,Vreal'*Gbus*Vreal + Vimag'*Gbus*Vimag)
    optimize!(model)
    return value.(Vn), value.(Pg), value.(Qg)
end

function LoadFlow(Ybus,Sn)
    n = length(Ybus[:,1])
    Ynn = Ybus[2:n,2:n]
    Yns = Ybus[2:n,1]
    Vn  = ones(n)*(1+0im)
    for k = 1:10 
        Vn[2:n] = Ynn\(conj.(Sn[2:n]./Vn[2:n])-Yns)
    end
    return Vn
end


function main()
    println("Reading data frames")
    Lines = CSV.read("Lines.csv",DataFrame)
    Nodes = CSV.read("Nodes.csv",DataFrame)
    Ybus = CalculateYbus(Lines)

    println("Calculating OPF using Wirtinger Linearization")
    Vr, Pg, Qg = LinearOPF(Lines,Nodes)
    Nodes.Pgen = round.(Pg,digits=4)
    Nodes.Qgen = round.(Qg,digits=4)
    Nodes.Vlin = abs.(Vr)
    Nodes.Alin  = angle.(Vr)*180/pi

    println("Evaluating power flow")
    Sn = Pg+Qg*1im-Nodes.Pload-Nodes.Qload*1im
    Vn = LoadFlow(Ybus,Sn)
    Nodes.Vmag = abs.(Vn)
    Nodes.Ang  = angle.(Vn)*180/pi

return Nodes
end

Nodes = main()
println(Nodes)
