num_nodes = 9;
num_lines = 9;
vec_lines = Array{NamedTuple,1}(undef, num_lines);
vec_lines[1] = (n1=1, n2=4, r=0.0000, x=0.0576, b=0.0000);
vec_lines[2] = (n1=2, n2=7, r=0.0000, x=0.0625, b=0.0000);
vec_lines[3] = (n1=3, n2=9, r=0.0000, x=0.0586, b=0.0000);
vec_lines[4] = (n1=4, n2=5, r=0.0100, x=0.0850, b=0.0880);
vec_lines[5] = (n1=4, n2=6, r=0.0170, x=0.0920, b=0.0790);
vec_lines[6] = (n1=5, n2=7, r=0.0320, x=0.1610, b=0.1530);
vec_lines[7] = (n1=6, n2=9, r=0.0390, x=0.1700, b=0.1790);
vec_lines[8] = (n1=7, n2=8, r=0.0085, x=0.0720, b=0.0745);
vec_lines[9] = (n2=8, n1=9, r=0.0119, x=0.1008, b=0.1045);

function Base.show(io::IO, vec_d::Vector{NamedTuple})
    printstyled(io,"Vector of NamedTuples of size ",length(vec_d),"\n";color=:green)
    w = vec_d[1]
    nd = length(w)
    q = keys(w)
    for k = 1:nd
        print(io,q[k],'\t')
    end
    print("\n")
    for d in vec_d
        for k in q
            if typeof(d[k])==Float64
               r = round(d[k],digits=4)
            else
               r = d[k]
            end
            print(io, "$r\t")
        end
        print("\n")
    end
end

 println(vec_lines)