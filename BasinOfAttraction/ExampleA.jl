# Estimación de la basija de atracción para el flujo de carga en una red DC simple

using LinearAlgebra
using Plots
using PolynomialRoots

vs = 1.0
pn  = [-0.5 -0.8]'
r12 = 0.12
r13 = 0.15
r23 = 0.13
Y = zeros(3,3)
Y[1,1] = 1/r12+1/r13
Y[1,2] =-1/r12
Y[1,3] =-1/r13
Y[2,1] =-1/r12
Y[2,2] = 1/r12+1/r23
Y[2,3] = -1/r23
Y[3,1] = -1/r13
Y[3,2] = -1/r23
Y[3,3] = 1/r13+1/r23

Zn = inv(Y[2:3,2:3])
Ym = Y[2:3,1]
T(v) = Zn*(pn./v-Ym*vs)
# flujo de carga
vf = zeros(2,10)
vf[:,1] = [0.4 0.4]'
for k = 1:9
    vf[:,k+1] = T(vf[:,k])
end
vff = vf[:,10]

vn = ones(2,1)
α1 = norm(Zn*pn,Inf)
α2 = norm(vn-Zn*(pn./vn-Ym*vs),Inf)

pol = [-α2;1+2*α2-α1;-2-α2;1]
r = sort(real(roots(pol)))   
rmin = r[1]
rmax = r[2]
println("vmin: ",minimum(vff)," vmin(est): ",1-rmin)
println("vmax: ",maximum(vff)," vmax(est): ",1+rmin)

theme(:dark)
plot([1-rmax,1-rmax,1+rmax,1+rmax,1-rmax],
      [1-rmax,1+rmax,1+rmax,1-rmax,1-rmax],
      color=:skyblue)
plot!([1-rmin,1-rmin,1+rmin,1+rmin,1-rmin],
      [1-rmin,1+rmin,1+rmin,1-rmin,1-rmin],
      color=:orange)
plot!([1],[1],marker=:circle,color=:orange)
anim = @animate for k = 1:10
    plot!(vf[1,1:k],vf[2,1:k],marker=:circle, 
          color=:green, xlims = (1-rmax-0.1,1+rmax+0.1), 
          ylims = (1-rmax-0.1,1+rmax+0.1), legend=false)
end 
gif(anim,fps=5)
