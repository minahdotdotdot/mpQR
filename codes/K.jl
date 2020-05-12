#mpTSQR2, change rows
include("mpTSQR.jl")
include("genmat.jl")
using DelimitedFiles, Printf

name = "K"
l = Float16; h=Float32; d=Float64;
ms = 2 .^(9:11); n = 100
trials = 10;
c = 10.
berr = zeros(d, length(ms), trials);
ferr = zeros(d, length(ms), trials);
L = 2;
writedlm("../txtfiles/"*name*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["mpBQR: forward error"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,m) in enumerate(ms)
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