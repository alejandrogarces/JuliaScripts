# Convergencia del flujo de carga
using XLSX 
using CSV
using DataFrames 
using Plots
using LinearAlgebra
using ProgressBars
using PolynomialRoots
using SparseArrays

# Tomar los datos de la hoja de excel
function LoadFeeder(file::String)
    DataLines = DataFrame(XLSX.readtable(file,"lines")) 
    Codes = DataFrame(XLSX.readtable(file,"line_codes")) 
    General = DataFrame(XLSX.readtable(file,"general")) 
    Loads = DataFrame(XLSX.readtable(file,"loads"))
    Profiles = DataFrame(XLSX.readtable(file,"profiles"))
    Pbase = General[1,2]/3           # potencia nominal por fase
    Vbase = General[1,1]/sqrt(3)     # voltaje linea a neutro
    Zbase = Vbase*Vbase/Pbase    
    NumN = maximum([maximum(DataLines[:,1]),maximum(DataLines[:,2])])
    NumL = length(DataLines[:,1])
    NumT = length(Profiles[:,1])
    Ybus = complex(zeros(3*NumN,3*NumN))
    # datos de linea
    Zline = complex(zeros(3,3))  
    for k = 1:NumL
        n1 = DataLines[k,1]
        n2 = DataLines[k,2]
        len = DataLines[k,3]/1000
        cde = DataLines[k,4]
        r1 = Codes[cde,2]*len/Zbase
        x1 = Codes[cde,3]*len/Zbase
        r0 = Codes[cde,4]*len/Zbase
        x0 = Codes[cde,5]*len/Zbase
        zs = complex((2*r1+r0)/3 , (2*x1+x0)/3)
        zm = complex((r0-r1)/3 , (x0-x1)/3)    
        Zline[1,1] = zs
        Zline[1,2] = zm
        Zline[1,3] = zm
        Zline[2,1] = zm
        Zline[2,2] = zs
        Zline[2,3] = zm
        Zline[3,1] = zm
        Zline[3,2] = zm
        Zline[3,3] = zs        
        yL = inv(Zline)
        nt1 = [n1,n1+NumN,n1+2*NumN]
        nt2 = [n2,n2+NumN,n2+2*NumN]    
        Ybus[nt1,nt1] = Ybus[nt1,nt1] + yL
        Ybus[nt1,nt2] = Ybus[nt1,nt2] - yL
        Ybus[nt2,nt1] = Ybus[nt2,nt1] - yL
        Ybus[nt2,nt2] = Ybus[nt2,nt2] + yL   
    end
    Vs = [cis(0),cis(-2*pi/3),cis(2*pi/3)]*General[1,3]
    # cargas
    Snode = complex(zeros(3*NumN,NumT))
    for load in eachrow(Loads)
        n = load[1]+(load[2]-1)*NumN
        Powerfactor = load[3]
        pro = load[4]
        Snode[n,:] = -Profiles[:,pro]*(1+1im*sqrt(1/Powerfactor^2-1))/Pbase/1000        
    end
return Ybus, Vs, Snode, NumN
end

# Flujo de carga
function LoadFlowFP(Ynn,Sn,Vn,In)    
    ew  = zeros(10)
    for iter = 1:10
        V = Ynn\(conj.(Sn./Vn)-In)
        err = norm(V-Vn)
        Vn = V        
        ew[iter] = err
    end
    return Vn, ew
end


function Bacin(Ybus, Vs, Sn, NumN)
    v0 = sum(abs.(Vs))/3
    Nslack = [1,NumN+1,2*NumN+1]    
    L = abs.(Sn)
    M = findall(k->k>0,L)
    MS = vcat(Nslack,M)  
    R = setdiff(1:3*NumN,MS)
    m = length(MS)
    YK = Ybus[MS,MS] - Ybus[MS,R]*(Ybus[R,R]\Ybus[R,MS])
    Ynn = YK[4:m,4:m]
    Znn = inv(Ynn)
    Yns = YK[4:m,1:3]
    V = kron(Vs./abs.(Vs)*v0,ones(NumN,1))
    Vn = V[M]
    Is = Znn*Yns*Vs+Vn        
    In = conj.(Sn[M]./Vn)
    Vnn = Znn*In
    α1 = norm(Znn*conj.(Sn[M]),Inf)
    α2 = norm(Vnn+Is,Inf)
    p = [-α2*v0^2;v0^2+2*α2*v0-α1;-2*v0-α2;1]
    r = sort(real(roots(p)))    
    rmin = r[1]
    rmax = r[2]    
    return α1,α2,rmin,rmax
end


# Main function
function main()
println("Loading data set")
Ybus, Vs, Snode, NumN = LoadFeeder("FEEDER900.xlsx")
t = 566
Nslack = [1,NumN+1,2*NumN+1]
Nother = setdiff(1:3*NumN,Nslack)
Sn = Snode[Nother,t]
Yss = Ybus[Nslack,Nslack]
Yns = Ybus[Nother,Nslack]
Ynn = Ybus[Nother,Nother]
Ynn = sparse(Ynn)
In = Yns*Vs
α1,α2,rmin,rmax = Bacin(Ybus, Vs, Sn, NumN)
Vini = kron(Vs,ones(NumN-1,1))
_, err = LoadFlowFP(Ynn,Sn,(1-rmax)*Vini,In)
plt = plot(err,yaxis=:log,color=:gray, legend = false)
VectErr = zeros(10,100)
for k = 1:100
    nu = (1-rmax)+k/100*(2*rmax)
    _, err = LoadFlowFP(Ynn,Sn,nu*Vini,In)
    plt = plot!(err,yaxis=:log,color=:gray, legend = false)
    VectErr[:,k] = err
end
_, err = LoadFlowFP(Ynn,Sn,Vini,In)    
plt = plot!(err,yaxis=:log,color=:green, marker = 4, legend = false,
            xlabel="iteración", ylabel="error")
theme(:dark)
display(plt)
return VectErr
end


Verr = main()
