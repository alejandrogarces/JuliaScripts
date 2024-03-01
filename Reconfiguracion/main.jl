using CSV
using DataFrames
using LinearAlgebra
using Plots

# Calcula la matriz Ybus a partir de las estructuras de datos
function calcular_ybus(datos,num_nodos,num_lineas)
    Ybus = zeros(num_nodos,num_nodos)*1im
    for k = 1:num_lineas
        n1 = datos.from[k]
        n2 = datos.to[k]
        y = 1/(datos.r_pu[k]+datos.x_pu[k]*1im)*datos.on[k]
        Ybus[n1,n1] +=  y 
        Ybus[n1,n2] += -y
        Ybus[n2,n1] += -y
        Ybus[n2,n2] +=  y
    end
    return Ybus
end

# Calcula el flujo de carga
function calcular_flujo_carga(Ybus,S,N,num_nodos,slack)
    YNN = Ybus[N,N]
    YNS = Ybus[N,slack]
    s_nodal = -S[N]
    v = ones(num_nodos,1)*(1+0im)
    err = 100
    iter = 0
    while err > 1E-8
        v[N] = YNN \ (conj.(s_nodal./v[N])-YNS*v[slack])  
        err = norm(s_nodal-v[N].*conj.(YNN*v[N]+YNS*v[slack]))
        iter += 1
        if iter > 100
            println("No converge despues de 100 iteraciones")
            break
        end
    end
    return v,err,iter
end

# funcion principal
function main()
    printstyled("Cargando los datos...\n", color = :green)
        datos = DataFrame(CSV.File("data_cinvalar.csv"))
        println(datos)
        num_nodos = maximum([maximum(datos.from),maximum(datos.to)])
        num_lineas = length(datos.from)
        Sbus = zeros(num_nodos,1)*1im
        slack = [1,2,3] # varias subestaciones
        N = setdiff(1:num_nodos,slack)
        for k = 1:num_lineas
            n = datos.to[k]
            Sbus[n] += (datos.p_MW[k]+datos.q_MW[k]*1im)/100;
        end
    printstyled("Calculando la Ybus...\n", color = :green)
        Ybus = calcular_ybus(datos,num_nodos,num_lineas)
        
    printstyled("Calculando flujo de carga...\n", color = :green)
        v,err,iter = calcular_flujo_carga(Ybus,Sbus,N,num_nodos,slack)    
        results = DataFrame(v_pu = abs.(v[:]), angle_deg = angle.(v[:])*180/pi)   
        println("Converge despues de ",iter," iteraciones")
        println(results)
    return 0
end

main()


