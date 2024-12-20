# Example of a port-controlled Hamiltonian System:
# Single Machine Connected to an Infinite Bus 

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

# Simulation with input = 0
results = Simulate(PCHS=SMIB,nt=5000,xini=[0.1;0.1],method="MidPoint")
PlotResults(PCHS=SMIB,data=results)
sleep(1)

# Simulation with input = 5
results = Simulate(PCHS=SMIB,nt=5000,xini=[0.1;0.1],u=[5],method="MidPoint")
PlotResults(PCHS=SMIB,data=results)
sleep(1)

# Simulation with input as step function
u_input(x,t) = [3*Heaviside(t-4)+5]
results = Simulate(PCHS=SMIB,nt=5000,xini=[0.1;0.1],u=u_input,method="MidPoint")
PlotResults(PCHS=SMIB,data=results)
sleep(1)

# Speed governor
tau = 2.5
kappa = 1
Ht(x) = 1/(2*tau)*x[1]^2
dHt(x) = x[1]/(tau)
GOV = BuildPCHS(nx=1,nu=1, dt=1/60/10, 
R=1, H=Ht, dH=dHt,G=kappa,
name="Speed Governor")

# Feedback
S = SMIB*GOV
u_in(x,t) = [-2*x[1];3*Heaviside(t-1)]
println(S)
results = Simulate(PCHS=S,nt=5000,u=u_in,method="MidPoint")
PlotResults(PCHS=S,data=results)

r = AdmissibleEquilibrium(PCHS=S)