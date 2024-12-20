include("scr/JuPCHS.jl")
using .JuPCHS
Ms = 30
wB = 2*pi*60
pmax = 10
Js = [0 -wB; wB 0]
Gs = [1 ; 0]
Hs(x) = 1/(2*Ms)*x[1]^2-pmax/wB*cos(x[2])
dHs(x) = [x[1]/Ms; pmax/wB*sin(x[2])]

println("Single machine infinite bus")
SMIB = BuildPCHS(nx=2,nu=1, dt=1/60/10, 
                 J=Js, H=Hs, dH=dHs,G=Gs,
                 name="SingleMachine",
                 xlabel=["x_ω","x_δ"],
                 ulabel=["p_mec"])
println(SMIB)

xm,hm = Hmin(PCHS=SMIB,xini = [0;pi/3])
println(xm)

AdmissibleEquilibrium(PCHS=SMIB,x=xm)