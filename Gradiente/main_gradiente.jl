##################################################################
####### Diferentes versiones del método del gradiente para #######
####### formas cuadráticas:  por Alejandro Garcés Ruiz.    #######
##################################################################
using LinearAlgebra
using Plots

# Imprimir para pgfplots
function pgf(x,y)
    n = length(x)
    s = ""
    for k = 1:n
        s *= "("
        s *=  string(x[k])
        s *= ","
        s *= string(y[k])
        s *= ")"
    end
    return s
end

# Gradiente con paso fijo
function GradientePasoFijo(f,g,t,xini,niter)
    x = xini
    gr_f = zeros(niter)
    gr_g = zeros(niter)
    gr_x = zeros((length(xini),niter))
    for k = 1:niter
        gr_f[k] = f(x)
        gr_g[k] = norm(g(x))
        gr_x[:,k] = x
        x += -t*g(x)
    end
    return gr_f,gr_g,gr_x
end


# Gradiente con paso de Polyak
function GradientePolyak(f,g,fopt,xini,niter)
    x = xini
    gr_f = zeros(niter)
    gr_g = zeros(niter)
    gr_x = zeros((length(xini),niter))
    for k = 1:niter
        gr_f[k] = f(x)
        gr_g[k] = norm(g(x))
        gr_x[:,k] = x
        t = (f(x)-fopt)/(norm(g(x)))^2
        x += -t*g(x)
    end
    return gr_f,gr_g,gr_x
end

function LineSearch(f, x, grad; t=1.0, tau=0.9)
    tp = t
    while f(x - tp * grad) > f(x)
        tp *= tau*tp
    end
    return tp
end

function GradienteLineSearch(f,g,t,xini,niter)    
    x = xini
    grad = g(x)
    grad_norm = norm(grad)
    gr_f = zeros(niter)
    gr_g = zeros(niter)
    gr_x = zeros((length(xini),niter))
    for k = 1:niter
        gr_f[k] = f(x)
        gr_g[k] = grad_norm
        gr_x[:,k] = x
        tp = LineSearch(f, x, grad, t=t)
        x += -tp*grad 
        grad = g(x)
        grad_norm = norm(grad)
    end
    return gr_f,gr_g,gr_x
end


# Gradiente con Backsteping
function GradientePolyak(f,g,fopt,xini,niter)
    x = xini
    gr_f = zeros(niter)
    gr_g = zeros(niter)
    gr_x = zeros((length(xini),niter))
    for k = 1:niter
        gr_f[k] = f(x)
        gr_g[k] = norm(g(x))
        gr_x[:,k] = x
        t = (f(x)-fopt)/(norm(g(x)))^2
        x += -t*g(x)
    end
    return gr_f,gr_g,gr_x
end


# Main function
Q = [5 -1;-1 2]
f(x) = 0.5*x'*Q*x
g(x) = Q*x
mu = norm(Q)
gamma = norm(Q)
t = 2*mu/(gamma^2)/2
xini = [1;1]
niter = 20
th = (1:50)*2*pi/100
U,Q2 = cholesky(Q)
Qi = inv(Q2)
y = [cos.(th),sin.(th)]
global plt_iter = plot([0 0],[0 0],label="",ylabel="y",xlabel="x")
for k = 1:5
    local yy = y*k/f(xini)
    local xx = Qi*yy
    global plt_iter = plot(plt_iter,xx[1],xx[2],label="",linecolor=:gray)
end
xini = ones(2)
str_paso = "Fijo t_opt"
println("Metodo del gradiente con paso ",str_paso)
gr_f,gr_g,gr_x = GradientePasoFijo(f,g,t,xini,niter)
plt_conv = plot(gr_g,yaxis=:log, marker = :circle, minorgrid=true,label=str_paso,ylabel="Grad",legend=:bottomleft)
plt_fobj = plot(gr_f,yaxis=:log, marker = :circle, minorgrid=true,label=str_paso,ylabel="f(x)",legend=:bottomleft)
plt_iter = plot(plt_iter,gr_x[1,:],gr_x[2,:],marker = :circle,label=str_paso)
#println(pgf(gr_x[1,:],gr_x[2,:]))

str_paso = "Fijo 1.9t_opt"
println("Metodo del gradiente con paso ",str_paso)
gr_f,gr_g,gr_x = GradientePasoFijo(f,g,1.9*t,xini,niter)
plt_conv = plot(plt_conv,gr_g,yaxis=:log, marker = :circle, label=str_paso)
plt_fobj = plot(plt_fobj,gr_f,yaxis=:log, marker = :circle, label=str_paso)
plt_iter = plot(plt_iter,gr_x[1,:],gr_x[2,:],marker = :circle,label=str_paso)
#println(pgf(gr_x[1,:],gr_x[2,:]))

str_paso = "Polyak"
println("Metodo del gradiente con paso ",str_paso)
gr_f,gr_g,gr_x = GradientePolyak(f,g,0,xini,niter)
plt_conv = plot(plt_conv,gr_g,yaxis=:log, marker = :circle,label=str_paso)
plt_fobj = plot(plt_fobj,gr_f,yaxis=:log, marker = :circle,label=str_paso)
plt_iter = plot(plt_iter,gr_x[1,:],gr_x[2,:],marker = :circle,label=str_paso)

str_paso = "line search"
println("Metodo del gradiente ",str_paso)
gr_f,gr_g,gr_x = GradienteLineSearch(f,g,4*t,xini,niter)
plt_conv = plot(plt_conv,gr_g,yaxis=:log, marker = :circle,label=str_paso)
plt_fobj = plot(plt_fobj,gr_f,yaxis=:log, marker = :circle,label=str_paso)
plt_iter = plot(plt_iter,gr_x[1,:],gr_x[2,:],marker = :circle,label=str_paso)
println(pgf(1:niter,gr_g))

# Graficos totales
plot!(plt_conv,plt_fobj,plt_iter, layout=[2,1],size = (600,600))
