include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")

using DelimitedFiles
name = "Y"
l = Float16; h=Float32; d=Float64;
m = 4096; ns = 2 .^(8:11);
trials = 10;
c=1000.;
berr = zeros(d, length(ns), trials,8);
oerr = zeros(d, length(ns), trials,8);
r = 128; 
L = 2;
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"f.txt", ["Condition orthogonal"])
for t = 1 : trials
	@printf("\nTrial%d: ",t)
	for (i,n) in enumerate(ns)
       @printf("n=%d, ",n)
        A = genmat(m,n,c,l);
        Ah = Matrix{h}(A);
        Ad = Matrix{d}(A);
       @printf("%d\t",1)

        # Block high precision
        Q, R = bhh_QR(Ah, r);
        berr[i,t,1] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,1] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

       @printf("%d\t",2)
        # mpBQR2
        Q, R = bhh_QR(A, r);
        berr[i,t,2] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,2] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

       @printf("%d\t",3)
        # mpBQR3
        Q, R = mpbhh_QR(A, r);
        berr[i,t,3] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,3] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        if m/n<4
           if m/n<2
              L=0;
           else
              L=1
           end
        else
          L=2
        end
        # TSQR high precision
       @printf("%d\t",4)
        Q, R = par_TSQR(Ah, L);
        berr[i,t,4] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,4] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR2
       @printf("%d\t",5)
        Q, R = par_TSQR(A, L);
        berr[i,t,5] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,5] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR3
       @printf("%d\t",6)
        Q, R = mpTSQR(A, L);
        berr[i,t,6] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,6] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

       @printf("%d\t",7)
        # hhQR in high precision
        Q, R = hh_QR(Ah);
        berr[i,t,7] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,7] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

       @printf("%d\t",8)
        # hhQr in low precision
        Q, R = hh_QR(A);
        berr[i,t,8] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,8] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
        
    end
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:,t,1]'])
        writedlm(io, [berr[:,t,2]'])
        writedlm(io, [berr[:,t,3]'])
        writedlm(io, [berr[:,t,4]'])
        writedlm(io, [berr[:,t,5]'])
        writedlm(io, [berr[:,t,6]'])
        writedlm(io, [berr[:,t,7]'])
        writedlm(io, [oerr[:,t,8]'])
        writedlm(io, ['\n'])
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [oerr[:,t,1]'])
        writedlm(io, [oerr[:,t,2]'])
        writedlm(io, [oerr[:,t,3]'])
        writedlm(io, [oerr[:,t,4]'])
        writedlm(io, [oerr[:,t,5]'])
        writedlm(io, [oerr[:,t,6]'])
        writedlm(io, [oerr[:,t,7]'])
        writedlm(io, [oerr[:,t,8]'])
        writedlm(io, ['\n'])
    end
end
