# Ejemplo de estabilidad transitoria usando el modelo clasico
# Basado en Power System Control and Stability 2nd Edition 
# Paul M. Anderson, A. A. Fouad

"""
    kron_reduction(Y::Matrix{ComplexF64},n:Int64)

TBW
"""
function kron_reduction(Y::Matrix{ComplexF64},n::Int64)
    n_t = length(Y[:,1]);  # tamano total de la matriz
    n_elim = 1:(n_t-n);    # nodos a ser eliminados
    n_rema = (n_t-n+1):n_t;  # nodos que permanecen
    y_kron = Y[n_rema,n_rema] - Y[n_rema,n_elim]*inv(Y[n_elim,n_elim])*Y[n_elim,n_rema]; 
    println(y_kron);
    return y_kron;
end


function Base.show(io::IO, vec_d::Vector{NamedTuple})
    printstyled(io,"Vector of NamedTuples of size ",length(vec_d),"\n";color=:green)
    w = vec_d[1]
    nd = length(w)
    q = keys(w)
    for k = 1:nd
        print(io,q[k],'\t')
    end
    print("\n")
    for d in vec_d
        nd = length(d)
        for k = 1:nd
            if typeof(d[k])==Float64
               r = round(d[k],digits=4)
            else
               r = d[k]
            end
            print(io, "$r\t")
        end
        print("\n")
    end
end

num_nodes = 9;
num_lines = 9;
vec_lines = Array{NamedTuple,1}(undef, num_lines);
vec_lines[1] = (n1=1, n2=4, r=0.0000, x=0.0576, b=0.0000);
vec_lines[2] = (n1=2, n2=7, r=0.0000, x=0.0625, b=0.0000);
vec_lines[3] = (n1=3, n2=9, r=0.0000, x=0.0586, b=0.0000);
vec_lines[4] = (n1=4, n2=5, r=0.0100, x=0.0850, b=0.0880);
vec_lines[5] = (n1=4, n2=6, r=0.0170, x=0.0920, b=0.0790);
vec_lines[6] = (n1=5, n2=7, r=0.0320, x=0.1610, b=0.1530);
vec_lines[7] = (n1=6, n2=9, r=0.0390, x=0.1700, b=0.1790);

vec_lines[8] = (n1=7, n2=8, r=0.0085, x=0.0720, b=0.0745);
vec_lines[9] = (n1=8, n2=9, r=0.0119, x=0.1008, b=0.1045);
println("lineas:",vec_lines)

num_loads = 3
vec_loads = Array{NamedTuple,1}(undef,num_loads);
vec_loads[1] = (n=5, p=1.25, q=0.50);
vec_loads[2] = (n=6, p=0.90, q=0.30);
vec_loads[3] = (n=8, p=1.00, q=0.35);
println("Cargas:",vec_loads);

num_gen = 3;
vec_gen = Array{NamedTuple,1}(undef,num_gen);
vec_gen[1] = (n=1, x=0.0608, H=23.640);
vec_gen[2] = (n=2, x=0.1198, H=6.8000);
vec_gen[3] = (n=3, x=0.1813, H=3.0100);
println("Generadores:",vec_gen);

# Flujo de carga
vn = [1.0400;1.0250;1.0250;1.0258;0.9956;1.0127;1.0258;1.0159;1.0324];
an = [0;9.2800;4.6648;-2.2168;-3.9888;-3.6874;3.7197;0.7275;1.9667];                  
v_node = vn.*exp.(an*pi/180*1im);

# Calculo de la ybus
y_bus = complex(zeros(num_nodes,num_nodes));
for line in vec_lines
    local y_line, b_cap
    y_line = 1/(line.r + line.x*1im);
    b_cap  = line.b*1im; 
    y_bus[line.n1,line.n1] = y_bus[line.n1,line.n1] + y_line + b_cap;
    y_bus[line.n1,line.n2] = y_bus[line.n1,line.n2] - y_line;
    y_bus[line.n2,line.n1] = y_bus[line.n2,line.n1] - y_line;
    y_bus[line.n2,line.n2] = y_bus[line.n2,line.n2] + y_line + b_cap;
end

i_node = y_bus*v_node
s_node = v_node.*conj.(i_node)
# agregar las cargas en la ybus
printstyled("Modelo de las cargas\n",color=:green);
for load in vec_loads
    z_load = v_node[load.n]*conj(v_node[load.n]/(load.p+load.q*1im));
    y_bus[load.n,load.n] = y_bus[load.n,load.n] + 1/z_load;
    println('z',load.n," : ",z_load)
end
# voltaje interno de los generadores
# ybus aumentada prefalla
n_aumentada = num_nodes+num_gen;
y_aumentada = complex(zeros(n_aumentada,n_aumentada));
N = 1:num_nodes
y_aumentada[N,N] = y_bus;
v_gen = zeros(num_gen,1);
a_gen = zeros(num_gen,1);

printstyled("Voltaje interno de los generadores\n",color=:green)
for k in 1:num_gen
    gen = vec_gen[k]
    vg = v_node[gen.n] + i_node[gen.n]*gen.x*1im;
    v_gen[k] = abs(vg);
    a_gen[k] = angle(vg);
    n1 = gen.n;
    n2 = gen.n + num_nodes;
    y_aumentada[n1,n1] = y_aumentada[n1,n1] + 1/(gen.x*1im);
    y_aumentada[n1,n2] = y_aumentada[n1,n2] - 1/(gen.x*1im);
    y_aumentada[n2,n1] = y_aumentada[n2,n1] - 1/(gen.x*1im);
    y_aumentada[n2,n2] = y_aumentada[n2,n2] + 1/(gen.x*1im);
    println('v',gen.n," : ", abs(vg)," < " ,angle(vg)*180/pi, "p : ", real(s_node[gen.n]));
end

# reduccion de nodos de Kron
printstyled("Matriz de prefalla\n",color=:green)
y_prefalla = kron_reduction(y_aumentada,num_gen);

# Matriz Falla
printstyled("Matriz en falla\n",color=:green)
n_t = setdiff(1:n_aumentada,7); # nodos que no estan en falla
YA = y_aumentada[n_t,n_t];
y_falla = kron_reduction(YA,num_gen)

# Matriz de postfalla
printstyled("Matriz en post-falla\n",color=:green)
YB = y_aumentada;
line = vec_lines[6]
y_line = 1/(line.r + line.x*1im);
b_cap  = line.b*1im; 
YB[line.n1,line.n1] = YB[line.n1,line.n1] - y_line - b_cap;
YB[line.n1,line.n2] = YB[line.n1,line.n2] + y_line;
YB[line.n2,line.n1] = YB[line.n2,line.n1] + y_line;
YB[line.n2,line.n2] = YB[line.n2,line.n2] - y_line - b_cap;
y_post = kron_reduction(YB,num_gen)

# valor de M
w_base = 2*pi*60
M = [2*gen.H for gen in vec_gen]
p_mec = [real(s_node[gen.n]) for gen in vec_gen]
# solucion de la ecuacion diferencial
t_falla = 0.1;
t_apertura = 0.183;
t_fin = 2.0;
num_steps = 60;
dt = t_fin/num_steps;
gr_delta = zeros(num_steps,num_gen);
gr_omega = zeros(num_steps,num_gen);
gr_time = (1:num_steps) * dt;

w  = ones(num_gen);
for k = 1:num_steps
    if gr_time[k] < t_falla
        Ybus = y_prefalla;
    elseif gr_time[k] < t_apertura
        Ybus = y_falla;
    else
        Ybus = y_post;
    end
    vg = v_gen.*exp.(a_gen*1im);
    ig = Ybus*vg;
    sg = vg.*conj(ig);
    pg = real(sg);
    for i = 1:num_gen
        w[i] = w[i] + dt/M[i]*(p_mec[i]-pg[i]);
        a_gen[i] = a_gen[i] + dt*(w[i]-1)*w_base;
        gr_delta[k,i] = a_gen[i]*180/pi;
        gr_omega[k,i] = w[i];
    end
end

using Plots
p1 = plot(gr_time,gr_delta)
p2 = plot(gr_time,gr_omega)
plot(p1,p2,layout=(2,1))

# sacar los datos en pgfplots
#i = 1
#Mt = sum(M)
#for k = 1:num_steps
#    t = round(gr_time[k],digits=4)
#    ci = sum(gr_delta[k,:].*M)/Mt
#    th = round(gr_delta[k,i]-ci,digits=4)   
#    print("($t,$th)")
#end