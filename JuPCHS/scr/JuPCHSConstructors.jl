"""
# Construct port hamiltonian system and check coherence #
# function BuildPCHS(;nx=1,nu=1,name="build_PHS",
              J=nothing,R=nothing,
              dH=nothing,H=nothing,G=nothing,u=nothing,
              xlabel=nothing,ulabel=nothing,ylabel=nothing, dt=nothing) #
"""

# matrix to matrixfunction
function matrix_to_matrixfuntion(M=nothing)
    if (typeof(M)== Matrix{Int64})||(typeof(M)== Matrix{Float64})            
        Mm = x ->  M*1.0            
    end   
    if (typeof(M)==Vector{Int64})||(typeof(M)==Vector{Float64})
        nx = length(M)
        Mw = zeros(nx,1)
        Mw[:,1] = M
        Mm = x -> Mw
    end 
    if (typeof(M)==Int64)||(typeof(M)==Float64)
        Mw = zeros(1,1)
        Mw[1,1] = M
        Mm = x -> Mw
    end            
    if length(methods(M))==1            
        Mm = x -> M(x)
    end    
    return Mm
end


# Check properties of matrix J
function check_matrix_J(J=nothing,nx=1)
    if J===nothing        
        print_warning("J is a zero matrix")        
        Jm = x->  zeros((nx,nx))        
    else        
        Jm = matrix_to_matrixfuntion(J)
    end
    # check the correct size of J
    x = rand(nx)
    Jx = Jm(x)
    m = size(Jx)    
    if m[1] != m[2]
        print_warning("J is not square") 
    end
    if m[1] != nx
        print_warning("The size of J is different from nx")
    end
    # check if it is antysymmetric    
    k = 0
    while k < MAXITER
        k += 1
        xp = rand(nx)
        Jx = Jm(xp)
        Jt = transpose(Jx)            
        if norm(Jx+Jt)>EPSILON
           print_warning("J is not skew-symmetric")            
           k = 11
        end    
    end        
    return Jm
end

# Check properties of matrix R
function check_matrix_R(R=nothing,nx=1)
    if R===nothing
        print_warning("R is a zero matrix")
        Rm = x->  zeros((nx,nx))
    else        
        Rm = matrix_to_matrixfuntion(R)
    end
    x = rand(nx)
    Rx = Rm(x)
    m = size(Rx)    
    if m[1] != m[2]
        print_warning("R is not square") 
    end
    if m[1] != nx
        print_warning("The size of R is different from nx")
    end
    # check if it semidefinite
    k = 0
    while k < MAXITER
        k += 1
        xp = rand(nx)
        Rx = Rm(xp)        
        L = minimum(real(eigvals(Rx)))        
        if L < 0
           print_warning("R is not semidefinite") 
           k = 11
        end
    end    
    return Rm
end

# check properties of matrix G
function check_matrix_G(G=nothing,nx=1,nu=1)
    if G===nothing
        print_warning("G is a zero matrix")
        Gm = x->  zeros((nx,nu))
    else        
        Gm = matrix_to_matrixfuntion(G)
    end
    x = rand(nx)    
    Gx = Gm(x)
    m = size(Gx)       
    if (m[1] != nx)||(m[2] != nu)
        print_warning("The size of G is different from (nx,nu)")
    end
    return Gm
end

# check consistency of label x
function check_label_x(xlabel=nothing,nx=1)
    if xlabel===nothing
       lx = Vector{String}(undef,nx)
       for k = 1:nx
           lx[k] = "x"*string(k)
       end
    else
       if length(xlabel) != nx
          print_warning("Length of xlabel is not consistent")
       end 
       lx = xlabel
    end
    return lx
end

# check consistency of label u
function check_label_u(ulabel=nothing,nu=1)
    if ulabel===nothing
       lu = Vector{String}(undef,nu)
       for k = 1:nu
           lu[k] = "u"*string(k)
       end
    else
       if length(ulabel) != nu
          print_warning("Length of ulabel is not consistent")
       end 
       lu = ulabel
    end
    return lu
end

# check consistency of label y
function check_label_y(ylabel=nothing,nu=1)
    if ylabel===nothing
       ly = Vector{String}(undef,nu)
       for k = 1:nu
           ly[k] = "y"*string(k)
       end
    else
       if length(ylabel) != nu
          print_warning("Length of ylabel is not consistent")
       end 
       ly = ylabel
    end
    return ly
end


function line_integral(dH::Function,nx::Integer,dt::Float64,points::Integer)
    xa = zeros(nx)    
    f = 0.0
    for k=1:points
        t = dt*k
        xb = ones(nx)*t
        ga = sum(dH(xa))
        gb = sum(dH(xb))
        gm = sum(dH(0.5*xa+0.5*xb))
        f = f+dt/6*(ga+4*gm+gb)
        xa = xb
    end
    return f
end

function check_hamiltonian(H=nothing,dH=nothing,nx=1,dt=0.001)  
    if dH === nothing
        dHm = x -> zeros(nx)
    else
        dHm = x -> dH(x)
    end
    if H === nothing
        Hm = x -> 0.0
    elseif (typeof(H)== Matrix{Int64})||(typeof(H)== Matrix{Float64})
        Hm = x -> 1/2*x'*H*x
        dHm = x -> H*x
    else
        Hm = x -> H(x)
    end
    # check consistency of H and dH by a linear integral
    lint = line_integral(dHm,nx,dt,100)    
    xa = zeros(nx)
    xb = ones(nx)*100*dt
    if abs(Hm(xb)-Hm(xa)-lint)>EPSILON
       print_warning("∇H(x) is not consistent with H(x)")
    end
    return Hm, dHm
end

# main constructor for the type build_PHS
function BuildPCHS(;nx=1,nu=1,name="build_PHS",
              J=nothing,R=nothing,
              dH=nothing,H=nothing,G=nothing,u=nothing,
              xlabel=nothing,ulabel=nothing,ylabel=nothing, dt=nothing)
    println("Building a port-controlled Hamiltonian-System")          
    if (nx<=0)||(nx===nothing)
       nxp = 1 
       print_warning("nx must be greater than 0. Set to nx=1")
    else
       nxp = Int(nx)        
    end   
    if (nu<=0)||(nu===nothing)
        nup = 1 
        print_warning("nu must be greater than 0. Set to nx=1")
    else
        nup = Int(nu)
    end              
    if dt===nothing    
        dt = EPSILON
    end
    Jp = check_matrix_J(J,nxp)    
    Rp = check_matrix_R(R,nxp)
    Gp = check_matrix_G(G,nxp,nup)    
    Hp, dHp = check_hamiltonian(H,dH,nx,dt)
    label_xp = check_label_x(xlabel,nx)
    label_up = check_label_u(ulabel,nu)          
    label_yp = check_label_y(ylabel,nu)          
    p = typePCHS(Jp,Rp,dHp,Hp,Gp,nxp,nup,dt,label_xp,label_up,label_yp,name)    
    return p
end
