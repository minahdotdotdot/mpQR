#mpTSQR2, change levels
include("TSQR.jl")
include("genmat.jl")
using DelimitedFiles, Printf

name = "L"
l = Float16; h=Float32; d=Float64;
m = 4000; n = 100
trials = 10;
c = 10.
Ls = 1:5;

berr = zeros(d, length(Ls), trials);
ferr = zeros(d, length(Ls), trials);
writedlm("../txtfiles/"*name*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["mpBQR: forward error"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,L) in enumerate(Ls)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ad = Matrix{d}(A);

        # Block mixed precision
        Q, R = par_TSQR(A, L);
        berr[i,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)/norm(Ad)
        ferr[i,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)
    end

    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:, t]'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [ferr[:, t]'])
    end
end
