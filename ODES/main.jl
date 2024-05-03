using LinearAlgebra
using Plots

struct pHsystem
    name::String
    J::Function
    R::Function
    dH::Function
    H::Function
    G::Function
    u::Function
    nx :: Int
    nu :: Int 
end

function Base.show(io::IO, p::pHsystem)
    printstyled(io,"Port Hamiltonian System\n", color=:blue)
    println(io,"   .num_n:\t $(p.nx) \t (number of state variables)")
    println(io,"   .num_l:\t $(p.nu) \t (number of inputs)")    
end   

function ExplicitEuler(S,xini,dt,nt)
    tode = (1:nt)*dt
    xode = zeros(S.nx,nt)
    uode = zeros(S.nu,nt)
    yode = zeros(S.nu,nt)
    Hode = zeros(nt)
    x = xini
    xode[:,1] = x
    uode[:,1] = [S.u(x)]
    yode[:,1] = [transpose(S.G(x))*S.dH(x)]
    Hode[1] = S.H(x)
    for k = 2:nt
        f = (S.J(x)-S.R(x))*S.dH(x)+S.G(x)*S.u(x)
        x = x + dt*f
        xode[:,k] = x
        uode[:,k] = [S.u(x)]
        yode[:,k] = [transpose(S.G(x))*S.dH(x)]
        Hode[k] = S.H(x)    
    end
    return tode,xode,uode,yode,Hode
end


function ImplicitEuler(S,xini,dt,nt)
    tode = (1:nt)*dt
    xode = zeros(S.nx,nt)
    uode = zeros(S.nu,nt)
    yode = zeros(S.nu,nt)
    Hode = zeros(nt)
    x = xini
    xode[:,1] = x
    uode[:,1] = [S.u(x)]
    yode[:,1] = [transpose(S.G(x))*S.dH(x)]
    Hode[1] = S.H(x)
    for k = 2:nt
        xa = x
        f = (S.J(xa)-S.R(xa))*S.dH(xa)+S.G(xa)*S.u(xa)
        xn = x + dt*f
        err = norm(xn-xa)
        xa = xn
        m = 0
        while (err>1E-8)&(m<50)
            f = (S.J(xa)-S.R(xa))*S.dH(xa)+S.G(xa)*S.u(xa)
            xn = x + dt*f
            err = norm(xn-xa)
            xa = xn
            m += 1
        end
        if m >= 50
            println("No converge el punto fijo")
        end 
        x = xn
        xode[:,k] = x
        uode[:,k] = [S.u(x)]
        yode[:,k] = [transpose(S.G(x))*S.dH(x)]
        Hode[k] = S.H(x)    
    end
    return tode,xode,uode,yode,Hode
end


function ExplicitRungeKutta(S,xini,dt,nt)
    tode = (1:nt)*dt
    xode = zeros(S.nx,nt)
    uode = zeros(S.nu,nt)
    yode = zeros(S.nu,nt)
    Hode = zeros(nt)
    x = xini
    xode[:,1] = x
    uode[:,1] = [S.u(x)]
    yode[:,1] = [transpose(S.G(x))*S.dH(x)]
    Hode[1] = S.H(x)
    F(x) = (S.J(x)-S.R(x))*S.dH(x)+S.G(x)*S.u(x)
    for k = 2:nt
        f1 = F(x)
        f2 = F(x+dt*f1/2)
        f3 = F(x+dt*f2/2)
        f4 = F(x+dt*f3)
        x = x + dt*(f1+2*f2+2*f3+f4)/6
        xode[:,k] = x
        uode[:,k] = [S.u(x)]
        yode[:,k] = [transpose(S.G(x))*S.dH(x)]
        Hode[k] = S.H(x)    
    end
    return tode,xode,uode,yode,Hode
end


function ImplicitRungeKutta(S,xini,dt,nt)
    tode = (1:nt)*dt
    xode = zeros(S.nx,nt)
    uode = zeros(S.nu,nt)
    yode = zeros(S.nu,nt)
    Hode = zeros(nt)
    x = xini
    xode[:,1] = x
    uode[:,1] = [S.u(x)]
    yode[:,1] = [transpose(S.G(x))*S.dH(x)]
    Hode[1] = S.H(x)
    F(x) = (S.J(x)-S.R(x))*S.dH(x)+S.G(x)*S.u(x)
    for k = 2:nt
        xa = x
        f1 = F(xa)
        f2 = F(xa+dt*f1/2)
        f3 = F(xa+dt*f2/2)
        f4 = F(xa+dt*f3)
        xn = x + dt*(f1+2*f2+2*f3+f4)/6
        err = norm(xn-xa)
        xa = xn
        m = 0
        while (err>1E-8)&(m<50)
            f1 = F(xa)
            f2 = F(xa+dt*f1/2)
            f3 = F(xa+dt*f2/2)
            f4 = F(xa+dt*f3)
            xn = x + dt*(f1+2*f2+2*f3+f4)/6            
            err = norm(xn-xa)
            xa = xn
            m += 1
        end
        if m >= 50
            println("No converge el punto fijo")
        end 
        x = xn
        xode[:,k] = x
        uode[:,k] = [S.u(x)]
        yode[:,k] = [transpose(S.G(x))*S.dH(x)]
        Hode[k] = S.H(x)    
    end
    return tode,xode,uode,yode,Hode
end


function main()
xini = [0.9;0.8]
tend = 1
dt = 1/60/4
nt = 800
wbase = 2*pi*60
M = 30
xi = 0
Pm = 5
Pmax = 10
J(x) = [0 -1; 1 0]
R(x) = diagm([xi/wbase,0])
H(x) = wbase/(2*M)*x[1]^2 + -Pmax*cos(x[2])
dH(x) = [wbase*x[1]/M; Pmax*sin(x[2])]
f(x)  = (J-R)*dH(x) + [1;0]*Pm
G(x) = [1;0]
u(x) = Pm
SMIB = pHsystem("Single Machine Infinite Bus System",J,R,dH,H,G,u,2,1)

theme(:dark)

tode,xode,uode,yode,Hode = ExplicitEuler(SMIB,xini,dt,nt)
p1 = plot(tode,xode[1,:],label="x_1 expEuler")
p2 = plot(tode,xode[2,:],label="x_2 expEuler")
p3 = plot(tode,uode[1,:],label="u expEuler")
p4 = plot(tode,Hode,label="H expEuler")

tode,xode,uode,yode,Hode = ImplicitEuler(SMIB,xini,dt,nt)
yu = yode.*uode
p1 = plot(p1,tode,xode[1,:],label="x_1 impEuler")
p2 = plot(p2,tode,xode[2,:],label="x_2 impEuler")
p3 = plot(p3,tode,uode[1,:],label="u impEuler")
p4 = plot(p4,tode,Hode,label="H impEuler")

tode,xode,uode,yode,Hode = ExplicitRungeKutta(SMIB,xini,dt,nt)
p1 = plot(p1,tode,xode[1,:],label="x_1 expRK")
p2 = plot(p2,tode,xode[2,:],label="x_2 expRK")
p3 = plot(p3,tode,uode[1,:],label="u expRK")
p4 = plot(p4,tode,Hode,label="H expRK")

tode,xode,uode,yode,Hode = ImplicitRungeKutta(SMIB,xini,dt,nt)
p1 = plot(p1,tode,xode[1,:],label="x_1 impRK")
p2 = plot(p2,tode,xode[2,:],label="x_2 impRK")
p3 = plot(p3,tode,uode[1,:],label="u impRK")
p4 = plot(p4,tode,Hode,label="H impRK")

plt = plot(p1,p2,p3,p4,layout=(2,2))
display(plt)
return SMIB
end


S = main()


