"""
# Numerical methods for ordinary differential equations #
# function simulate(;PCHS=nothing,xini=nothing,nt=nothing,method=nothing) #
# function explicit_euler(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64) #  
# function implicit_euler(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64) #
# function implicit_midle_point(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64) #
# function explicit_runge_kutta(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64) #
# function implicit_runge_kutta(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64) #
"""

# Explicit Euler: x+ = x + dt*F(x)
function explicit_euler(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64)    
    res = DataFrame()    
    res.t = (1:nt)*PCHS.dt
    xode = zeros(PCHS.nx,nt)
    uode = zeros(PCHS.nu,nt)
    yode = zeros(PCHS.nu,nt)    
    res.H = zeros(nt)
    x = xini
    ut = u(x,0)    
    xode[:,1] = x    
    uode[:,1] = ut
    yode[:,1] = transpose(PCHS.G(x))*PCHS.dH(x)
    res.H[1] = PCHS.H(x)    
    for k = 2:nt
        ut = u(x,(k-1)*PCHS.dt)
        f = (PCHS.J(x)-PCHS.R(x))*PCHS.dH(x)+PCHS.G(x)*ut
        x = x + PCHS.dt*f
        xode[:,k] = x
        uode[:,k] = ut
        yode[:,k] = transpose(PCHS.G(x))*PCHS.dH(x)
        res.H[k] = PCHS.H(x)    
    end   
    for k = 1:PCHS.nx 
        res[!,PCHS.xlabel[k]] = xode[k,:]
    end
    for k = 1:PCHS.nu 
        res[!,PCHS.ulabel[k]] = uode[k,:]
    end
    for k = 1:PCHS.nu
        res[!,PCHS.ylabel[k]] = yode[k,:]
    end
    return res
end

# Implicit Euler
function implicit_euler(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64)   
    res = DataFrame()    
    res.t = (1:nt)*PCHS.dt
    xode = zeros(PCHS.nx,nt)
    uode = zeros(PCHS.nu,nt)
    yode = zeros(PCHS.nu,nt)
    res.H = zeros(nt)    
    x = xini
    ut = u(x,0)
    xode[:,1] = x
    uode[:,1] = ut
    yode[:,1] = transpose(PCHS.G(x))*PCHS.dH(x)
    res.H[1] = PCHS.H(x)    
    for k = 2:nt
        xa = x
        ut = u(x,(k-1)*PCHS.dt)
        f = (PCHS.J(xa)-PCHS.R(xa))*PCHS.dH(xa)+PCHS.G(xa)*ut
        xn = x + PCHS.dt*f
        err = norm(xn-xa)
        xa = xn
        m = 0
        while (err>EPSILON)&(m<MAXITER)
            ut = u(xn,(k-1)*PCHS.dt)
            f = (PCHS.J(xa)-PCHS.R(xa))*PCHS.dH(xa)+PCHS.G(xa)*ut
            xn = x + PCHS.dt*f
            err = norm(xn-xa)
            xa = xn
            m += 1
        end
        if m >= MAXITER
            print_warning("Fixed-point algorithm does not converge")
        end 
        x = xn
        xode[:,k] = x
        uode[:,k] = ut
        yode[:,k] = transpose(PCHS.G(x))*PCHS.dH(x)
        res.H[k] = PCHS.H(x)            
    end
    for k = 1:PCHS.nx 
        res[!,PCHS.xlabel[k]] = xode[k,:]
    end
    for k = 1:PCHS.nu 
        res[!,PCHS.ulabel[k]] = uode[k,:]
    end
    for k = 1:PCHS.nu
        res[!,PCHS.ylabel[k]] = yode[k,:]
    end
    return res
end

# Implicit middle point
function implicit_midle_point(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64)   
    res = DataFrame()    
    res.t = (1:nt)*PCHS.dt
    xode = zeros(PCHS.nx,nt)
    uode = zeros(PCHS.nu,nt)
    yode = zeros(PCHS.nu,nt)
    res.H = zeros(nt)    
    x = xini
    ut = u(x,0)
    xode[:,1] = x
    uode[:,1] = ut
    yode[:,1] = transpose(PCHS.G(x))*PCHS.dH(x)
    res.H[1] = PCHS.H(x)    
    for k = 2:nt
        xm = x
        ut = u(xm,(k-1)*PCHS.dt)
        f = (PCHS.J(xm)-PCHS.R(xm))*PCHS.dH(xm)+PCHS.G(xm)*ut
        xm = x + PCHS.dt*f/2
        err = norm(xm-x)        
        m = 0
        while (err>EPSILON)&(m<MAXITER)
            xa = xm
            ut = u(xm,(k-1)*PCHS.dt)
            f = (PCHS.J(xm)-PCHS.R(xm))*PCHS.dH(xm)+PCHS.G(xm)*ut
            xm = x + PCHS.dt*f/2
            err = norm(xm-xa)            
            m += 1
        end
        if m >= MAXITER
            print_warning("Fixed-point algorithm does not converge")
        end 
        x = 2*xm-x
        xode[:,k] = x
        uode[:,k] = ut
        yode[:,k] = transpose(PCHS.G(x))*PCHS.dH(x)
        res.H[k] = PCHS.H(x)            
    end
    for k = 1:PCHS.nx 
        res[!,PCHS.xlabel[k]] = xode[k,:]
    end
    for k = 1:PCHS.nu 
        res[!,PCHS.ulabel[k]] = uode[k,:]
    end
    for k = 1:PCHS.nu
        res[!,PCHS.ylabel[k]] = yode[k,:]
    end
    return res
end

# Explicit RK
function explicit_runge_kutta(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64)    
    res = DataFrame()    
    res.t = (1:nt)*PCHS.dt
    xode = zeros(PCHS.nx,nt)
    uode = zeros(PCHS.nu,nt)
    yode = zeros(PCHS.nu,nt)    
    res.H = zeros(nt)
    x = xini
    ut = u(x,0)
    xode[:,1] = x    
    uode[:,1] = ut
    yode[:,1] = transpose(PCHS.G(x))*PCHS.dH(x)
    res.H[1] = PCHS.H(x)    
    F(x,ur) = (PCHS.J(x)-PCHS.R(x))*PCHS.dH(x)+PCHS.G(x)*ur
    for k = 2:nt
        ut = u(x,(k-1)*PCHS.dt)
        f1 = F(x,ut)
        f2 = F(x+PCHS.dt*f1/2,ut)
        f3 = F(x+PCHS.dt*f2/2,ut)
        f4 = F(x+PCHS.dt*f3,ut)
        x = x + PCHS.dt*(f1+2*f2+2*f3+f4)/6        
        xode[:,k] = x
        uode[:,k] = ut
        yode[:,k] = transpose(PCHS.G(x))*PCHS.dH(x)
        res.H[k] = PCHS.H(x)    
    end   
    for k = 1:PCHS.nx 
        res[!,PCHS.xlabel[k]] = xode[k,:]
    end
    for k = 1:PCHS.nu 
        res[!,PCHS.ulabel[k]] = uode[k,:]
    end
    for k = 1:PCHS.nu
        res[!,PCHS.ylabel[k]] = yode[k,:]
    end
    return res
end

# Implicit RK
function implicit_runge_kutta(PCHS::typePCHS ,u::Function,xini::Vector,nt::Int64)
    res = DataFrame()    
    res.t = (1:nt)*PCHS.dt
    xode = zeros(PCHS.nx,nt)
    uode = zeros(PCHS.nu,nt)
    yode = zeros(PCHS.nu,nt)
    res.H = zeros(nt)    
    x = xini
    ut = u(x,0)
    xode[:,1] = x
    uode[:,1] = ut
    yode[:,1] = transpose(PCHS.G(x))*PCHS.dH(x)
    res.H[1] = PCHS.H(x)    
    F(x,ur) = (PCHS.J(x)-PCHS.R(x))*PCHS.dH(x)+PCHS.G(x)*ur
    for k = 2:nt
        xa = x
        ut = u(x,(k-1)*PCHS.dt)
        f1 = F(xa,ut)
        f2 = F(xa+PCHS.dt*f1/2,ut)
        f3 = F(xa+PCHS.dt*f2/2,ut)
        f4 = F(xa+PCHS.dt*f3,ut)
        xn = x + PCHS.dt*(f1+2*f2+2*f3+f4)/6
        err = norm(xn-xa)
        xa = xn
        m = 0
        while (err>EPSILON)&(m<MAXITER)
            ut = u(xn,(k-1)*PCHS.dt)
            f1 = F(xa,ut)
            f2 = F(xa+PCHS.dt*f1/2,ut)
            f3 = F(xa+PCHS.dt*f2/2,ut)
            f4 = F(xa+PCHS.dt*f3,ut)
            xn = x + PCHS.dt*(f1+2*f2+2*f3+f4)/6            
            err = norm(xn-xa)
            xa = xn
            m += 1
        end
        if m >= MAXITER
            println("Fixed-point algorithm does not converge")
        end 
        x = xn
        xode[:,k] = x
        uode[:,k] = ut
        yode[:,k] = transpose(PCHS.G(x))*PCHS.dH(x)
        res.H[k] = PCHS.H(x)        
    end
    for k = 1:PCHS.nx 
        res[!,PCHS.xlabel[k]] = xode[k,:]
    end
    for k = 1:PCHS.nu 
        res[!,PCHS.ulabel[k]] = uode[k,:]
    end
    for k = 1:PCHS.nu
        res[!,PCHS.ylabel[k]] = yode[k,:]
    end
    return res
end


# Main function for simulation
function Simulate(;PCHS=nothing,u=nothing,xini=nothing,nt=nothing,method=nothing)
    res = DataFrame()
    if PCHS===nothing
        print_warning("Empty PCHS")
        PCHS =  build_PHS()
    end
    if xini===nothing
        xini = zeros(PCHS.nx) 
    end
    if nt === nothing 
       nt = MAXITER
    end    
    if u === nothing
       um = (x,t) -> zeros(PCHS.nu)
    else
        um = deepcopy(u)
        if (typeof(u)==Vector{Int64})||(typeof(u)==Vector{Float64})
            um = (x,t) -> u
        end 
        if (typeof(u)==Int64)||(typeof(u)==Float64)
            Mw = zeros(1,1)
            Mw[1,1] = u*1.0
            um = (x,t) -> Mw
        end  
    end
    if method in [nothing,"Euler","ForwardEuler"]
        res = explicit_euler(PCHS,um,xini,nt)    
    end
    if method in ["ImplicitEuler","BackwardEuler","iEuler","iE"]
        res = implicit_euler(PCHS,um,xini,nt)
    end            
    if method in ["MidPoint","ImplicitMidPoint","iMP"]
        res = implicit_midle_point(PCHS,um,xini,nt)
    end            
    if method in ["RungeKutta","RK"]
        res = explicit_runge_kutta(PCHS,um,xini,nt)
    end            
    if method in ["ImplicitRungeKutta","iRK"]
        res = implicit_runge_kutta(PCHS,um,xini,nt)
    end
    return res
end