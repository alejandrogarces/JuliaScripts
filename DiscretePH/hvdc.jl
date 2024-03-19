## Sistema HVDC modelo simplificado 
using LinearAlgebra
using Plots

# Parametros
l = 0.1
c = 1.1
r = 0.1
p = 1.2
v_bar = 1.0
# Sistema hamiltoniano controlado por puertos
J= [0 -1; 1 0]

phi(v) = 1/(1+exp(-100*(v-1.0)))
R(x) = [r 0; 0 phi(x[2]/c)]
H(x) = x[1]^2/(2*l) + x[2]^2/(2*c)
grad_H(x) = [x[1]/l; x[2]/c]
g = [1; 0]
u = 1.0
gamma = sqrt((1/l)^2+(1/c)^2)
lambda = p/v_bar^2
zeta = sqrt(2+r^2+p^2/v_bar^4)
tau = 1/(gamma*zeta)

function Backward_Euler_Method(x_ini,np)
    gr_x = zeros(2,np)
    gr_e = zeros(1,np)
    gr_t = zeros(1,np)
    x = x_ini
    for k = 1:np
        z = x
        er = 10.0
        while er > 1E-11
            T = x + tau*((J-R(z))*grad_H(z)+g*u)
            er = norm(T-z)
            z = T
        end
        y = tau*transpose(grad_H(z))*g
        gr_t[1,k] = tau
        x = z;
        gr_x[:,k] = x
        gr_e[k] = er
    end
return gr_x, gr_t, gr_e   
end

function Proposed_Euler_Method(x_ini,np)
    gr_x = zeros(2,np)
    gr_e = zeros(1,np)
    gr_t = zeros(1,np)
    x = x_ini
    for k = 1:np
        z = x
        er = 10.0
        while er > 1E-11
            T = x + tau*((J-R(x))*grad_H(z)+g*u)
            er = norm(T-z)
            z = T
        end
        y = tau*transpose(grad_H(z))*g
        gr_t[1,k] = tau
        x = z;
        gr_x[:,k] = x
        gr_e[k] = er
    end
return gr_x, gr_t, gr_e  
end


function Forward_Euler_Method(x_ini,np)
    gr_x = zeros(2,np)
    x = x_ini
    for k = 1:np
        z = x + tau*((J-R(x))*grad_H(x)+g*u)
        y = tau*transpose(grad_H(z))*g
        x = z;
        gr_x[:,k] = x
    end
return gr_x   
end

function export_to_tikz(v,np)
    s = " "
    for k = 1:np
        s = s * "("
        s = s * string(k)
        s = s * ","
        s = s * string(v[k])
        s = s * ")"
    end
    return s
end

function main()
np = 30
x_ini = [0;1]
gr_x,_,_  = Backward_Euler_Method(x_ini,np)
gr_x2 = Forward_Euler_Method(x_ini,np)
gr_x3,_,_ = Proposed_Euler_Method(x_ini,np)

theme(:dark)
plt1 = plot([gr_x[1,:],gr_x2[1,:],gr_x3[1,:]],
            linetype=:steppre, 
            title="x_1",
            labels=["Backward" "Forward" "Proposed"],
            legend=:bottomright)
plt2 = plot([gr_x[2,:],gr_x2[2,:],gr_x3[2,:]],
            linetype=:steppre, 
            title="x_2",
            labels=["Backward" "Forward" "Proposed"],
            legend=:bottomright)

gr_H1 = [H(gr_x[:,k]) for k in 1:np]
gr_H2 = [H(gr_x2[:,k]) for k in 1:np]
gr_H3 = [H(gr_x3[:,k]) for k in 1:np]

plt3 = plot([gr_H1,gr_H2,gr_H3],
            linetype=:steppre, 
            title="Hamiltonian",
            labels=["Backward" "Forward" "Proposed"],
            legend=:bottomright)

plt = plot(plt1,plt2,plt3)
display(plt)
return "by Alejandro Garces"
end

main()
