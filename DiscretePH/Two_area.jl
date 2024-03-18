## Sistema de dos areas en forma Hamiltoniana
using LinearAlgebra
using Plots

# Parameters
M1 = 1
M2 = 2
ξ1 = 0.5
ξ2 = 0.8
Pmax = 8.3
Pm1 = 3.5
Pm2 = -Pm1
w_base = 2*pi*60

# pH structure
J = [0 0 -1; 0 0 1; 1 -1 0]
R = [ξ1 0 0; 0 ξ2 0; 0 0 0]
H(x) = x[1]^2/(2*M1) + x[2]^2/(2*M2) - Pmax*cos(x[3])
grad_H(x) = [x[1]/M1; x[2]/M2; Pmax*sin(x[3])]
g = [1 0; 0 1; 0 0]
u = [Pm1;Pm2] 
gamma = norm([1/M1;1/M2;Pmax])
tau = 1/(gamma*norm(J-R))*0.99999

function Backward_Euler_Method(x_ini,np)
    gr_x = zeros(3,np)
    gr_e = zeros(1,np)
    gr_s = zeros(2,np)
    x = x_ini
    for k = 1:np
        z = x
        er = 10.0
        while er > 1E-11
            T = x + tau*((J-R)*grad_H(z)+g*u)
            er = norm(T-z)
            z = T
        end
        y = tau*transpose(grad_H(z))*g
        gr_s[1,k] = y*u
        gr_s[2,k] = H(z)-H(x)
        x = z;
        gr_x[:,k] = x
        gr_e[k] = er
    end
return gr_x, gr_s, gr_e   
end


function Forward_Euler_Method(x_ini,np)
    gr_x = zeros(3,np)
    gr_s = zeros(2,np)
    x = x_ini
    for k = 1:np
        z = x + tau*((J-R)*grad_H(x)+g*u)
        y = tau*transpose(grad_H(z))*g
        gr_s[1,k] = y*u
        gr_s[2,k] = H(z)-H(x)
        x = z;
        gr_x[:,k] = x
    end
return gr_x, gr_s    
end
    
function main()
np = 100
x_ini = [0;0;0.7]
gr_x, gr_s, gr_e = Backward_Euler_Method(x_ini,np)
gr_x2, gr_s2 = Forward_Euler_Method(x_ini,np)

plt1 = plot([gr_x[1,:],gr_x2[1,:]],
            linetype=:steppre, 
            title="x_1", 
            labels=["Backward" "Forward"])
plt2 = plot([gr_x[2,:],gr_x2[2,:]],
            linetype=:steppre,
            title="x_2",
            labels=["Backward" "Forward"])
plt3 = plot([gr_x[3,:],gr_x2[3,:]],
            linetype=:steppre,
            title="x_3",
            labels=["Backward" "Forward"])
gr_H1 = [H(gr_x[:,k]) for k in 1:np]
gr_H2 = [H(gr_x2[:,k]) for k in 1:np]
plt4 = plot([gr_H1,gr_H2],
            linetype=:steppre, 
            title="Hamiltonian",
            labels=["Backward" "Forward"])
plt = plot(plt1,plt2,plt3,plt4)
theme(:dark)
display(plt)
return "by Alejandro Garces"
end

main()
