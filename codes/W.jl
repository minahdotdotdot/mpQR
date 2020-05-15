include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")

using DelimitedFiles
name = "W3"
l = Float16; h=Float32; d=Float64;
#c=1000.;
#m=4096;n=2048;
#rs = 2 .^(6:10); 
#berr = zeros(d, length(rs), trials,3);
#oerr = zeros(d, length(rs), trials,3);
m=2^13
n=m
rs = 2 .^ (12:-1:2)
trials = 1;
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
for t = 1 : trials
	#@printf("\n%d: ",t)
	for (i,r) in enumerate(rs)
        @printf("%d\t",i)
        open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,['\n'])
        end
        Ah = randn(h, m, n); Ah=Ah/norm(Ah)
        A = Matrix{l}(Ah);#genmat(m,n,c,l);
        #Ah = Matrix{h}(Ad);
        #Ad = Matrix{d}(A);

        # Block high precision
        Q, R = bhh_QR(Ah, r);
        b1 = norm(Matrix{h}(Q)*Matrix{h}(R)-Ah)
        f1 = opnorm(Matrix{h}(Q')*Matrix{h}(Q)-I)
        open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,[b1])
        end
        open("../txtfiles/"*name*"f.txt", "a") do io
            writedlm(io,[f1])
        end

        #=
        # mpBQR2
        Q, R = bhh_QR(A, r);
        b2 = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        f2 = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
        open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,[b2])
        end
        open("../txtfiles/"*name*"f.txt", "a") do io
            writedlm(io,[f2])
        end
        =#
print("mp3\t")
        # mpBQR3
        Q, R = mpbhh_QR(A, r);
        b3 = norm(Matrix{h}(Q)*Matrix{h}(R)-Ah)
        f3 = opnorm(Matrix{h}(Q')*Matrix{h}(Q)-I)
        open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,[b3])
        end
        open("../txtfiles/"*name*"f.txt", "a") do io
            writedlm(io,[f3])
        end
    end
    #=open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:,t,1]'])
        writedlm(io, [berr[:,t,2]'])
        writedlm(io, [berr[:,t,3]'])
        writedlm(io, ['\n'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [oerr[:,t,1]'])
        writedlm(io, [oerr[:,t,2]'])
        writedlm(io, [oerr[:,t,3]'])
        writedlm(io, ['\n'])
    end
    =#
end
