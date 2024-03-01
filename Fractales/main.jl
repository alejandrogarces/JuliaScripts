using LinearAlgebra
function h(x)
    theta_2 = x[1]
    theta_3 = x[2]
    theta_23 = theta_2-theta_3
    h1 = 0.4923 - 0.3978*cos(theta_2)-1.888*sin(theta_2)-0.2439*cos(theta_23)-2.0154*sin(theta_23)
    h2 = 3.9234 - 1.6390*cos(theta_3)-7.269*sin(theta_3)-0.2439*cos(theta_23)+2.0154*sin(theta_23)
    return [h1,h2]
end

function Jac(x)
    theta_2 = x[1]
    theta_3 = x[2]
    theta_23 = theta_2-theta_3
    J11 =  0.3978*sin(theta_2)-1.888*cos(theta_2)+0.2439*sin(theta_23)-2.0154*cos(theta_23) # dh1/dtheta_2
    J12 = -0.2439*sin(theta_23)+2.0154*cos(theta_23)
    J21 =  0.2439*sin(theta_23)+2.0154*cos(theta_23)
    J22 =  1.6390*sin(theta_3)-7.269*cos(theta_3)-0.2439*sin(theta_23)-2.0154*cos(theta_23)
    return [J11 J12;J21 J22]
end   

function NewtonRaphson(xini,h,Jac)
    x = xini
    hx = h(x)
    er = Inf
    k = 0
    while (er>1e-8)&(k<20)        
        Jx = Jac(x)
        x += - Jx\hx
        hx = h(x)
        er = norm(hx)
        k += 1      
    end
    return [er,k]   
end

xmax = pi
np = 2*314
p = range(-xmax,xmax,length=np)
frac_er = zeros(np,np)
frac_it = zeros(np,np)
for k in 1:np
    for m in 1:np
        xi = [p[k],p[m]]
        a,b = NewtonRaphson(xi,h,Jac)        
        frac_er[k,m] = a
        frac_it[k,m] = b        
    end
end
using Plots
plt = heatmap(frac_it)
display(plt)
