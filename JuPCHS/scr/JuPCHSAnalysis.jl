"""
Tools for analysis
# function Heaviside(t)
# function AdmissibleEquilibrium(;PCHS=nothing,u=nothing,x=nothing)
# function Hmin(;PCHS=nothing,xini=nothing)
"""

function Heaviside(t)
    return 0.5 * (sign(t) + 1)
end

function AdmissibleEquilibrium(;PCHS=nothing,u=nothing,x=nothing)
    if PCHS===nothing
        print_warning("Empty PCHS")
        PCHS =  build_PHS()
    end
    if u===nothing
        print_warning("Assuming input u=0")
        u = zeros(PCHS.nu)
    end
    if x===nothing
        print_warning("Assuming states x=0")
        x = zeros(PCHS.nx)
    end
    A = (PCHS.J(x)-PCHS.R(x))
    w = (A*PCHS.dH(x) + PCHS.G(x)*u)/norm(A)
    r = norm(w)
    rs = true
    if r > EPSILON
        print_warning("It is not an admissible equilibrium")
        println("|Δ| =",r, ", dx=", w)
        rs = false
    end
    return rs
end

function line_search(f, g, x, p; α_init=1.0, c1=1e-4, tau=0.5)
    α = α_init
    while f(x + α * p) > f(x) + c1 * α * (g(x)' * p)
        α *= tau
    end
    return α
end

function nonlinear_conjugate_gradient(f, g, x0; tol=EPSILON, max_iter=MAXITER)
    x = x0
    grad = g(x)
    p = -grad
    grad_norm = norm(grad)
    for k = 1:max_iter
        if grad_norm < tol
            break
        end
        α = line_search(f, g, x, p)
        x_new = x + α * p
        grad_new = g(x_new)
        β = max(0, (grad_new' * (grad_new - grad)) / (grad' * grad))
        p = -grad_new + β * p
        x = x_new
        grad = grad_new
        grad_norm = norm(grad)
    end
    return x, grad_norm
end



function Hmin(;PCHS=nothing,xini=nothing)
    if xini === nothing
        xini = zeros(PCHS.nx)
    end
    println("Optimizing by the conjugate gradient method")
    x_opt, grad_norm = nonlinear_conjugate_gradient(PCHS.H, PCHS.dH, xini)
    println("|∇H| =", grad_norm)
    hm = PCHS.H(x_opt)
    return x_opt, hm
end



