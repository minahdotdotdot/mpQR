#FIX n to 1000.
using MAT, DelimitedFiles
include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")


l=Float16; h=Float32; d=Float64
ms = exp10.(range(3, stop=log10(14000), length=10));
ns = 250 *ones(10);
ms = ceil.(Int, ms);
L = 2;
κs = 1e3*ones(length(ms));
dt = h;
name="size3_"
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
function addtofile!(name::String, b::AbstractFloat, o::AbstractFloat)
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io,[b])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io,[o])
    end
end

for i = 2 : 10#length(ms)
	varname="st3"*string(i,pad=2)
	file = matopen("../"*varname*".mat")
	A = Matrix{l}(read(file, varname))
	close(file)
	Ah = Matrix{h}(A)

    Q, R = hh_QR(A);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);
    
    n = size(A)[2];
    r = ceil(Int, n/4);
    # Block high precision
    Q, R = bhh_QR(Ah, r);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # mpBQR2
    Q, R = bhh_QR(A, r);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # mpBQR3
    Q, R = mpbhh_QR(A, r);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # TSQR high precision
    Q, R = par_TSQR(Ah, L);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # mpTSQR2
    Q, R = par_TSQR(A, L);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # mpTSQR3
    Q, R = mpTSQR(A, L);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);

    # hhQr
    Q, R = hh_QR(Ah);
    Q = Matrix{h}(Q);
    b = norm(Q*Matrix{h}(R)-Ah)
    o= opnorm(Q'*Q-I)
    addtofile!(name, b, o);
    
    addtofile!(name, -1., -1.);
    print("Done with"* string(ms[i])*"!\n")
end



function savetestmats!(ms::Vector{Int}, ns::Vector{Int}, κs,
	dt::DataType, name::String)
	for i = 1 : length(ms)
		varname=name*string(i,pad=2)
		file = matopen("../"*varname*".mat", "w")
		write(file, varname, genmat(ms[i],ns[i],κs[i],dt))
		close(file)
	end
end
