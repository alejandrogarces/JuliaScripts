using LinearAlgebra
using Plots

struct pHsystem
    name::String
    J::Matrix
    R::Matrix
    dH::Function
    H::Function
    g::Function
    u::Function 
end

f(x) =  [1-2*sin(x[2])-10*x[1]; 2*pi*60*x[1]]

# Metodo de Euler
function ode_Euler(f,x_ini,dt,steps)
    t = zeros(steps)
    n = length(x_ini)
    x = zeros(n,steps)
    h = zeros(steps)
    x[:,1] = x_ini
    h[1] = pi*60*x_ini[1]^2 - 2*cos(x_ini[2])
    for k = 1:steps-1
        x[:,k+1] = x[:,k] + f(x[:,k])*dt
        t[k+1] = t[k] + dt
        h[k+1] = pi*60*x[1,k+1]^2 - 2*cos(x[2,k+1])
    end
    return t,x,h
end

# Metodo de Runge-Kuta
function ode_RungeKuta(f,x_ini,dt,steps)
    t = zeros(steps)
    n = length(x_ini)
    x = zeros(n,steps)
    h = zeros(steps)
    x[:,1] = x_ini
    h[1] = pi*60*x_ini[1]^2 - 2*cos(x_ini[2])
    for k = 1:steps-1
        K1 = f(x[:,k])
        K2 = f(x[:,k]+0.5*dt*K1)
        K3 = f(x[:,k]+0.5*dt*K2)
        K4 = f(x[:,k]+dt*K3)
        x[:,k+1] = x[:,k] + dt*(K1+2*K2+2*K3+K4)/6
        t[k+1] = t[k] + dt
        h[k+1] = pi*60*x[1,k+1]^2 - 2*cos(x[2,k+1])
    end
    return t,x,h
end



steps = 300
dt = 0.01
x_ini = [0.0;0.0]
t_euler,x_euler,h_euler = ode_Euler(f,x_ini,dt,steps)
t_rk,x_rk,h_rk = ode_RungeKuta(f,x_ini,dt,steps)

p1 = plot(t_euler,[x_euler[1,:],x_rk[1,:]])
p2 = plot(t_euler,[x_euler[2,:],x_rk[2,:]])
p3 = plot(t_euler,[h_euler,h_rk])

plot(p1,p2,p3,layout=(3,1))


