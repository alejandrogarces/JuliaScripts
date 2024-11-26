# importa de cdf a dataframes
using DataFrames
using CSV

function read_cdf(file_name)
    println("Reading "*file_name)
    s_base = 0.0
    f = open(file_name,"r")        
    # linea 1
    s = readline(f)
    s_base = parse(Float64,s[31:36])
    println("|- S_BASE = ", s_base)
    println("|- Reading BUS DATA")
    # linea 2: no tiene informacion relevante
    s = readline(f)
    # linea 3 en adelante: BUS DATA
    s = readline(f)
    NUMBER = Int64[]
    NAME = String[]
    TYPE = Int64[]
    VPU = Float64[]
    ANG = Float64[]
    PLOAD = Float64[]
    QLOAD = Float64[]
    PGEN = Float64[]
    QGEN = Float64[]
    kVBASE = Float64[]
    QMIN = Float64[]
    QMAX = Float64[]
    GNOD = Float64[]
    BNOD = Float64[]    
    while s[1:4] != "-999"                
        push!(NUMBER,parse(Int64,s[1:4]))    
        push!(NAME,s[6:17])
        push!(TYPE,parse(Int64,s[25:26]))
        push!(VPU,parse(Float64,s[28:33]))
        push!(ANG,parse(Float64,s[34:40]))
        push!(PLOAD,parse(Float64,s[41:49]))
        push!(QLOAD,parse(Float64,s[50:59]))
        push!(PGEN,parse(Float64,s[60:67]))
        push!(QGEN,parse(Float64,s[68:75]))
        push!(kVBASE,parse(Float64,s[77:83]))
        push!(QMIN,parse(Float64,s[99:106]))
        push!(QMAX,parse(Float64,s[91:98]))
        push!(GNOD,parse(Float64,s[107:114]))
        push!(BNOD,parse(Float64,s[115:122]))
        s = readline(f)
    end
    nodes = DataFrame()    
    nodes.NUMBER = NUMBER
    nodes.NAME = NAME
    nodes.TYPE = TYPE
    nodes.VPU = VPU
    nodes.ANG = ANG
    nodes.PLOAD = PLOAD/s_base
    nodes.QLOAD = QLOAD/s_base
    nodes.PGEN = PGEN/s_base
    nodes.QGEN = QGEN/s_base
    nodes.KVBASE = kVBASE
    nodes.QMIN = QMIN/s_base
    nodes.QMAX = QMAX/s_base
    nodes.G = GNOD
    nodes.B = BNOD
    nodes.ID = [k for k = 1:nrow(nodes)]
    println("|- Number of nodes: ",nrow(nodes))
    # lineas
    s = readline(f)
    if s[1:6] == "BRANCH"
        println("|- Reading BRANCH DATA")
    end
    s = readline(f)
    FROM = Int64[]
    TO = Int64[]
    TIPO = Int64[]
    R = Float64[]
    X = Float64[]
    B = Float64[]
    LIM1 = Float64[]
    LIM2 = Float64[]
    LIM3 = Float64[]
    TAP = Float64[]
    TAPMAX = Float64[]
    TAPMIN = Float64[]
    STEP = Float64[]
    TANG = Float64[]
    while (s[1:4] != "-999") 
        push!(FROM,parse(Int64,s[1:4]))    
        push!(TO,parse(Int64,s[6:9]))    
        push!(TIPO,parse(Int64,s[19]))    
        push!(R,parse(Float64,s[20:29]))    
        push!(X,parse(Float64,s[30:40]))    
        push!(B,parse(Float64,s[41:50]))    
        push!(LIM1,parse(Float64,s[51:55]))    
        push!(LIM2,parse(Float64,s[57:61]))    
        push!(LIM3,parse(Float64,s[63:67]))    
        push!(TAP,parse(Float64,s[77:82]))    
        push!(TAPMAX,parse(Float64,s[98:104]))    
        push!(TAPMIN,parse(Float64,s[91:97]))    
        push!(STEP,parse(Float64,s[106:111]))    
        push!(TANG,parse(Float64,s[84:90]))    
        s = readline(f)
    end
    lines = DataFrame()
    lines.FROM = FROM
    lines.TO = TO
    lines.TYPE = TIPO
    lines.R = R
    lines.X = X
    lines.B = B
    lines.LIM1 = LIM1/s_base
    lines.LIM2 = LIM2/s_base
    lines.LIM3 = LIM3/s_base
    lines.TAP = TAP
    lines.TAPMIN = TAPMIN
    lines.TAPMAX = TAPMAX
    lines.STEP = STEP
    lines.TANG = TANG
    close(f)
    println("|- Number of lines: ",nrow(lines))    
    println("|- Renaming the nodes")
    for k = 1:nrow(lines)        
        lines.FROM[k] = nodes[nodes.NUMBER .== lines.FROM[k], "ID"][1]
        lines.TO[k] = nodes[nodes.NUMBER .== lines.TO[k], "ID"][1]        
    end
    println("|- Reference angle slack")
    th_slack = nodes[nodes.TYPE .== 3, "ANG"][1]
    nodes.ANG = nodes.ANG.-th_slack
    println("End")
    return nodes,lines
end    

file_name = "IEEE300.cdf"
nodes,lines = read_cdf(file_name)
CSV.write("/nodes.csv",nodes)
CSV.write("/lines.csv",lines)
