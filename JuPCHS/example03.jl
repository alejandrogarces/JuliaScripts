# Series an parallel interconection
include("scr/JuPCHS.jl")
using .JuPCHS
S1 = BuildPCHS(name="S1", nx = 5, ulabel=["r"], dt=0.1)
S2 = BuildPCHS(name="S2", nx = 3)
S3 = S1+S2
S4 = S1*S2
S5 = 2*S1