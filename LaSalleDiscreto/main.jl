using LinearAlgebra

function dif_ecuacion()
    np = 20
    a = 1/2
    b = 1/10
    x = [2.2,0]
    gr_x = zeros(2,np)
    gr_x[:,1] = x
    for k = 2:np
        d = x[1]^2+x[2]^2+b
        x = a*[(1-1/d)*x[2];(-1+1/d)*x[1]] 
        gr_x[1,k] = x[1]
        gr_x[2,k] = x[2] 
    end
    return gr_x
end

A =
gr_x = dif_ecuacion()
l = 0.233
s = 0.345
th = 0:0.1:2*pi+0.1
r = sqrt(l+s)
plot(r*sin.(th),r*cos.(th)) 
anim = @animate for k in 1:length(gr_x[1,:])-1
    plot!([gr_x[1,k],gr_x[1,k+1]],[gr_x[2,k],gr_x[2,k+1]],marker=:circle, color=:green, xlims = (-1,1.5), ylims = (-1,0.8), legend=false)
end
gif(anim, fps=5)
