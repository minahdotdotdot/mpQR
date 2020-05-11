#mpBQR3, change condition number
include("blockHQR.jl")
include("mpblockHQR.jl")
include("genmat.jl")
using DelimitedFiles, Printf

name = "H"
l = Float16; h=Float32; d=Float64;
m = 1024; ns = 2 .^(6:10)
trials = 10;
c = 10.
berr = zeros(d, length(ns), trials);
ferr = zeros(d, length(ns), trials);
r = 64;
writedlm("../txtfiles/"*name*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["mpBQR: forward error"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,n) in enumerate(ns)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ad = Matrix{d}(A);

        # Block mixed precision
        Q, R = mpbhh_QR(A, r);
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