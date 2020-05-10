using PyPlot, Printf, DelimitedFiles
include("blockHQR.jl")
include("mpblockHQR.jl")

name = "BQR-r-5k-"
mpname ="mpBQR-r-5k"
l = Float16; h=Float32; d=Float64;
m = 5000;
rs = [100, 250, 500, 1000, 2500];
trials = 10;
berr = zeros(d, length(rs), 3, trials);
ferr = zeros(d, length(rs), 3, trials);
n = 5000;
writedlm("../txtfiles/"*name*"b.txt", ["BQR: backward error"])
writedlm("../txtfiles/"*name*"f.txt", ["BQR: forward error"])
writedlm("../txtfiles/"*mpname*"b.txt", ["mpBQR: backward error"])
writedlm("../txtfiles/"*mpname*"f.txt", ["mpBQR: forward error"])

for (i,r) in enumerate(rs)
    #@printf("\n n=%d, trial: ",n)
    
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [string(r)])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [string(r)])
    end
    open("../txtfiles/"*mpname*"b.txt", "a") do io
        writedlm(io, [string(r)])
    end
    open("../txtfiles/"*mpname*"f.txt", "a") do io
        writedlm(io, [string(r)])
    end
    
    for t = 1 : trials
       @printf("%d\t",t)
        A = randn(l, m, n);
        Ah = Matrix{h}(A);
        Ad = Matrix{d}(A);

        # Inner product mixed precision
        #Q, R = bhh_QR(A, r);
        #berr[i,1,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)
        #ferr[i,1,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)

        # Block mixed precision
        Q, R = mpbhh_QR(A, r);
        berr[i,2,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)/norm(Ad)
        ferr[i,2,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)


        # High precision 
        Q, R = bhh_QR(Ah, r);
        berr[i,3,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)/norm(Ad)
        ferr[i,3,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)
    end

    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[i, 3, :]'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [ferr[i, 3, :]'])
    end
    open("../txtfiles/"*mpname*"b.txt", "a") do io
        writedlm(io, [berr[i, 2, :]'])
    end
    open("../txtfiles/"*mpname*"f.txt", "a") do io
        writedlm(io, [ferr[i, 2, :]'])
    end
end
