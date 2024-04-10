using LinearAlgebra
using CSV
using DataFrames
using Plots
using JuMP
using Ipopt
using ProgressBars

function calculate_gbus(Data,n)
    G = zeros(n,n)
    p = zeros(n,1)
    for line in eachrow(Data)
        g = 1/line.rkm
        G[line.from,line.from] +=  g
        G[line.from,line.to]   += -g
        G[line.to,line.from]   += -g
        G[line.to,line.to]     +=  g
        p[line.to] += line.pk
    end
    return p,G
end

function load_flow(G,p,n)
    N = 2:n
    Gnn = G[N,N]
    Gns = G[N,1]
    v = ones(n,1)
    for k = 1:10
        v[N] = Gnn\(p[N]./v[N]-Gns*1)
    end
    return v
end

function calculate_Z(ψ,G,u,n)
    Z = zeros(n,n)    
    for i = 1:n, m = 1:n, k = 1:n        
        Z[i,m] += ψ[k]*G[k,i]*G[k,m]*u[k]^2                
    end      
    return Z
end


function calculate_y(ϕ,ψ,G,u,v_bar,p_bar,n)
    y = zeros(n)
    for i = 1:n
        y[i] = ϕ[i]*v_bar[i]
        for k = 1:n
            y[i] += ψ[k]*G[k,i]*u[k]*p_bar[k]
        end
    end
    return y
end

function calculate_T(ϕ,ψ,G,u,v_bar,p_bar,n)
    Z = calculate_Z(ψ,G,u,n)
    y = calculate_y(ϕ,ψ,G,u,v_bar,p_bar,n)
    un = (diagm(ϕ)+Z)\y  # Map T        
    return un
end

function Algorithm1(v_bar,p_bar,ϕ,ψ,G,n)
    np = 4
    u = ones(n)
    gr_s = zeros(np)    
    for k = 1:np
        un = calculate_T(ϕ,ψ,G,u,v_bar,p_bar,n)      
        gr_s[k] = norm(u-un)        
        u = un
    end
    display(plot(gr_s,yaxis=:log))    
    return u
end

function optimization_model(n,ϕ,ψ,p_bar,v_bar,G)
    model = Model(Ipopt.Optimizer)
    @variable(model, vs[1:n])
    @variable(model, ps[1:n])
    @objective(model,Min,0.5*sum((ϕ[k]*(vs[k]-v_bar[k])^2+ψ[k]*(ps[k]-p_bar[k])^2) for k = 1:n))
    for k = 1:n
        @constraint(model, ps[k] == vs[k]*sum(G[k,m]*vs[m] for m=1:n))    
        @constraint(model, abs(vs[k]-1) <= 0.2)
    end
    set_silent(model)
    optimize!(model)
    println(termination_status(model))
    vss = zeros(n)
    for k = 1:n
        vss[k] = value(vs[k])
    end    
    return vss    
end


function Algorithm2(ϕ,ψ,G,v_bar,p_bar,n)
    L = 0
    un = ones(n)
    println("Estimating Lipschits constant")
    for k in ProgressBar(1:1E5)
        u = un + 0.2*randn(n)
        v = un + 0.2*randn(n)
        tu = calculate_T(ϕ,ψ,G,u,v_bar,p_bar,n)
        tv = calculate_T(ϕ,ψ,G,v,v_bar,p_bar,n)
        d  = norm(tu-tv)/norm(u-v)
        L  = maximum([d,L]) 
    end
    return L
end


# main function
function main()
    printstyled("A Fixed-point Algorithm for the Non-linear State Estimation in DC Microgrids\n"; color = :blue)
    Data = DataFrame(CSV.File("microgrid.csv"))    
    n = maximum([maximum(Data.from),maximum(Data.to)])
    pn,G = calculate_gbus(Data,n)
    v = load_flow(G,pn,n)
    p = v.*(G*v)
    # random noise
    v_bar = v + randn(n)*0.001
    v_bar[1] = 1
    p_bar = p + randn(n)*0.01
    ϕ = 100*ones(n)
    ψ =   1*ones(n)
    # malfunction devices    
    n1 = Integer(round(rand()*(n-1)+1))
    v_bar[n1] = v[n1] + randn()*0.1
    ϕ[n1] = 1
    n2 = Integer(round(rand()*(n-1)+1))
    v_bar[n2] = v[n2] + randn()*0.1
    ϕ[n2] = 1
    # state estimation
    @time vs = optimization_model(n,ϕ,ψ,p_bar,v_bar,G)   
    @time u = Algorithm1(v_bar,p_bar,ϕ,ψ,G,n)    
    # printing results
    results = DataFrame()
    results.v_bareasured = zeros(n)
    results.v_fixpoint = zeros(n)
    results.v_optimum  = zeros(n)
    results.v_real = zeros(n)    
    for k = 1:n
        results.v_bareasured[k] = v_bar[k]
        results.v_fixpoint[k] = u[k]
        results.v_optimum[k]  = vs[k]
        results.v_real[k] = v[k]
    end
    results.err_fix = results.v_fixpoint-results.v_real
    results.err_opt = results.v_optimum-results.v_real    
    println(results)    
    L = Algorithm2(ϕ,ψ,G,v_bar,p_bar,n)
    println("L = ",L)    
   return "alejandro.garces@utp.edu.co"
end

main()
