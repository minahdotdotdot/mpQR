include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")

using DelimitedFiles
name = "sBQR3"
l = Float16; h=Float32; d=Float64;
#c=1000.;
#m=4096;n=2048;
#rs = 2 .^(6:10); 
#berr = zeros(d, length(rs), trials,3);
#oerr = zeros(d, length(rs), trials,3);
#m=2^13
#n=m
#rs = 2 .^ (12:-1:2)
m=2^11;
n=2^9;
rs = 2 .^(1:9);
trials = 10;
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
for t = 1 : trials
    Ah = genmat(m,n,1000.0,h);
    A = Matrix{l}(Ah);
    Ah = Matrix{h}(Ah)
	@printf("\n%d: ",t)
	for (i,r) in enumerate(rs)
        @printf("%d\t",i)
        #= open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,['\n'])
        end =#
        #Ah = randn(h, m, n); Ah=Ah/norm(Ah)
        #A = Matrix{l}(Ah);#genmat(m,n,c,l);
        #=open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,[string(cond(Ah))*"\n"])
        end=#
        
        
        # Block high precision
        Q, R = bhh_QR(Ah, r);
        b1 = norm(Matrix{h}(Q)*Matrix{h}(R)-Ah)
        f1 = opnorm(Matrix{h}(Q')*Matrix{h}(Q)-I)
        #=
        open("../txtfiles/"*name*"b.txt", "a") do io
            writedlm(io,[b1])
        end
        open("../txtfiles/"*name*"f.txt", "a") do io
            writedlm(io,[f1])
        end
        
        # mpBQR2
        Q, R = bhh_QR(A, r); Qh=Matrix{h}(Q);
        b2 = norm(Qh*Matrix{h}(R)-Ah)
        f2 = opnorm(Qh'*Qh-I)
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
            writedlm(io,[b1 b3])
        end
        open("../txtfiles/"*name*"f.txt", "a") do io
            writedlm(io,[f1 f3])
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
