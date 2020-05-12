#mpBQR3, change condition number
include("blockHQR.jl")
include("mpblockHQR.jl")
include("genmat.jl")
using DelimitedFiles, Printf

name = "C2"
mpname="mp"*name
l = Float16; h=Float32; d=Float64;
m = 4096; ns = 4096
trials = 10;
c = 100.
rs = 2 .^(6:10)
berr = zeros(d, length(rs), trials,2);
ferr = zeros(d, length(rs), trials,2);

writedlm("../txtfiles/"*name*"b.txt", ["BQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["BQR: forward error"])
writedlm("../txtfiles/"*mpname*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*mpname*"f.txt", ["mpBQR: forward error"])
for t = 1 : trials
    @printf("\n%d: ",t)
    for (i,r) in enumerate(rs)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ah = Matrix{h}(A);
        Ad = Matrix{d}(A);

        # Block mixed precision
        Q, R = bhh_QR(Ah, r);
        berr[i,t,1] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)/norm(Ad)
        ferr[i,t,1] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)

        # Block mixed precision
        Q, R = mpbhh_QR(A, r);
        berr[i,t,2] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)/norm(Ad)
        ferr[i,t,2] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)
    end
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:, t,1]'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [ferr[:, t,1]'])
    end

    open("../txtfiles/"*mpname*"b.txt", "a") do io
        writedlm(io, [berr[:, t,2]'])
    end
    open("../txtfiles/"*mpname*"f.txt", "a") do io
        writedlm(io, [ferr[:, t],2'])
    end
end