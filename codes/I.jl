#mpBQR3, change condition number
include("blockHQR.jl")
include("mpblockHQR.jl")
include("genmat.jl")
using DelimitedFiles, Printf

name = "I"
l = Float16; h=Float32; d=Float64;
m = 1024; n = 1024;
trials = 10;
c = 10.
rs = 2 .^(5:9);


berr = zeros(d, length(rs), trials);
ferr = zeros(d, length(rs), trials);

writedlm("../txtfiles/"*name*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["mpBQR: forward error"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,r) in enumerate(rs)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ad = Matrix{d}(A);

        # Block mixed precision
        Q, R = bhh_QR(A, r);
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