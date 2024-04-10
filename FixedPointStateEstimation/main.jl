using LinearAlgebra
using CSV
using DataFrames
using Plots
using JuMP
using Ipopt

function calculate_gbus(Data,num_l,num_n)
    G = zeros(num_n,num_n)
    p = zeros(num_n,1)
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

function load_flow(G,p,num_n)
    N = 2:num_n
    Gnn = G[N,N]
    Gns = G[N,1]
    v = ones(num_n,1)
    for k = 1:10
        v[N] = Gnn\(p[N]./v[N]-Gns*1)
    end
    return v
end

function calculate_Z(phi,G,u,num_n)
    Z = zeros(num_n,num_n)
    for i = 1:num_n
        for m = 1:num_n
            for k = 1:num_n
                Z[i,m] += phi[k]*G[k,i]*G[k,m]*u[k]^2
            end
        end
    end
    return Z
end

function calculate_y(phi,psi,G,u,v_m,p_m,num_n)
    y = zeros(num_n)
    for i = 1:num_n
        y[i] = phi[i]*v_m[i]
        for k = 1:num_n
            y[i] += psi[k]*G[k,i]*u[k]*p_m[k]
        end
    end
    return y
end

function fixed_point(v_m,p_m,phi,psi,G,num_n)
    np = 3
    u = ones(num_n)
    gr_s = zeros(np)
    for k = 1:np
        Z = calculate_Z(phi,G,u,num_n)
        y = calculate_y(phi,psi,G,u,v_m,p_m,num_n)
        u = (diagm(phi)+Z)\y  # Map T        
        pu = u.*(G*u)
        gr_s[k] = 0.5*sum(phi.*(u-v_m).^2) + 0.5*sum(psi.*(pu-p_m).^2)        
    end
    #display(plot(gr_s,yaxis=:log))
    println(gr_s)
    return u
end

function optimization_model(num_n,phi,psi,p_m,v_m,G)
    model = Model(Ipopt.Optimizer)
    @variable(model, vs[1:num_n])
    @variable(model, ps[1:num_n])
    @objective(model,Min,sum((phi[k]*(vs[k]-v_m[k])^2+psi[k]*(ps[k]-p_m[k])^2) for k = 1:num_n))
    for k = 1:num_n
        @constraint(model, ps[k] == vs[k]*sum(G*vs))
    end
    set_silent(model)
    optimize!(model)
    println(termination_status(model))
    return vs
end

# main function
#function main()
    Data = DataFrame(CSV.File("microgrid.csv"))
    num_l = nrow(Data)
    num_n = maximum([maximum(Data.from),maximum(Data.to)])
    p,G = calculate_gbus(Data,num_l,num_n)
    v = load_flow(G,p,num_n)
    p = v.*(G*v)
    # random noise
    v_m = v + randn(num_n)*0.001
    v_m[1] = 1
    p_m = p + randn(num_n)*0.01
    phi = 100*ones(num_n)
    psi =   1*ones(num_n)
    # state estimation
    u = fixed_point(v_m,p_m,phi,psi,G,num_n)
    vs = optimization_model(num_n,phi,psi,p_m,v_m,G)
    # printing results
    results = DataFrame()
    results.v_measured = zeros(num_n)
    results.v_fixpoint = zeros(num_n)
    results.v_optimum  = zeros(num_n)
    results.v_real = zeros(num_n)
    for k = 1:num_n
        results.v_measured[k] = v_m[k]
        results.v_fixpoint[k] = u[k]
        results.v_optimum[k]  = value(vs[k])
        results.v_real[k] = v[k]
    end
    println(results)
#   return u
#end

#main()