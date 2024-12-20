
# Example of a port-controlled Hamiltonian System:
# RLC circuit

# Bulding the system
using LinearAlgebra
include("scr/JuPCHS.jl")
using .JuPCHS
R = diagm([1,5,1/50,1/2])
Q = diagm([1/1.3E-3,1/1.5E-3,1/200E-6,1/300E-6])
J = [0 0 -1 0; 0 0 1 -1; 1 -1 0 0; 0 1 0 0]
G = [1;0;0;0]
CIRCUIT = BuildPCHS(nx=4,nu=1,H=Q,J=J,R=R,G=G,dt=1E-5)
println(CIRCUIT)

results = Simulate(PCHS=CIRCUIT, u=[20])
PlotResults(PCHS=CIRCUIT,data=results,xscale=Q)