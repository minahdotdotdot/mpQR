using MAT
ms = exp10.(range(3, stop=log10(45000), length=14));
ns = floor.(Int,ms/4);
ms = ceil.(Int, ms);
L = 2;
κs = 1e3*ones(length(ms));
dt = h;
name="size"
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
function addtofile!(name::String, b::Float64, o::Float64)
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io,[b])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io,[o])
    end
end

for i = 1 : 9#length(ms)
	varname=name*string(i,pad=2)
	file = matopen("../"*varname*".mat")
	A = Matrix{l}(read(file, varname))
	close(file)
	Ah = Matrix{h}(A)
	Ad = Matrix{d}(A)

    Q, R = hh_QR(A);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o= opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
    addtofile!(name, b, o);
    
    n = size(A)[2];
    r = ceil(Int, n/4)
    # Block high precision
    Q, R = bhh_QR(Ah, r);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
    addtofile!(name, b, o);

    # mpBQR2
    Q, R = bhh_QR(A, r);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
	addtofile!(name, b, o);

    # mpBQR3
    Q, R = mpbhh_QR(A, r);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
	addtofile!(name, b, o);

    # TSQR high precision
    Q, R = par_TSQR(Ah, L);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
	addtofile!(name, b, o);

    # mpTSQR2
    Q, R = par_TSQR(A, L);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
	addtofile!(name, b, o);

    # mpTSQR3
    Q, R = mpTSQR(A, L);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
	addtofile!(name, b, o);

    # hhQr
    Q, R = hh_QR(Ah);
    b = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
    o = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
    addtofile!(name, b, o);
    
    addtofile!(name, -1., -1.);
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