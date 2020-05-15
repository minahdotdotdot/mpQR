include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")

using DelimitedFiles
name = "V"
l = Float16; h=Float32; d=Float64;
m = 4096; n=128;
trials = 10;
c=1000.;
Ls = 1:5;
berr = zeros(d, length(Ls), trials,5);
oerr = zeros(d, length(Ls), trials,5);
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,L) in enumerate(Ls)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ah = Matrix{h}(A);
        Ad = Matrix{d}(A);
        # TSQR high precision
        Q, R = par_TSQR(Ah, L);
        berr[i,t,1] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,1] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR2
        Q, R = par_TSQR(A, L);
        berr[i,t,2] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,2] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR3
        Q, R = mpTSQR(A, L);
        berr[i,t,3] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,3] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
        
        # hhQR in high precision
        Q, R = hh_QR(Ah);
        berr[i,t,4] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,4] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # hhQr in low precision
        Q, R = hh_QR(A);
        berr[i,t,5] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,5] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)    
        
    end
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:,t,1]'])
        writedlm(io, [berr[:,t,2]'])
        writedlm(io, [berr[:,t,3]'])
        writedlm(io, [berr[:,t,4]'])
        writedlm(io, [berr[:,t,5]'])
        writedlm(io, ['\n'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [oerr[:,t,1]'])
        writedlm(io, [oerr[:,t,2]'])
        writedlm(io, [oerr[:,t,3]'])
        writedlm(io, [berr[:,t,4]'])
        writedlm(io, [berr[:,t,5]'])
        writedlm(io, ['\n'])
    end
end
