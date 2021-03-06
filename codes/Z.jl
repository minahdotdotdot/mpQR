include("genmat.jl")
include("blockHQR.jl")
include("mpblockHQR.jl")
include("mpTSQR.jl")
include("TSQR.jl")

using DelimitedFiles
name = "Z2"
l = Float16; h=Float32; d=Float64;
m = 4096; n = 1024;
trials = 10;
cs = exp10.(range(0, stop=4, length=9));
berr = zeros(d, length(cs), trials)#,7);
oerr = zeros(d, length(cs), trials)#,7);
r = 256; 
L = 2;
writedlm("../txtfiles/"*name*"b.txt", ["Condition backward"])
writedlm("../txtfiles/"*name*"fo.txt", ["Condition orthogonal"])
for t = 1 : trials
	@printf("\n%d: ",t)
	for (i,c) in enumerate(cs)
       @printf("%d\t",i)
        A = genmat(m,n,c,l);
        Ah = Matrix{h}(A);
        Ad = Matrix{d}(A);
        # hhQR low precision
        Q, R = hh_QR(A);
        berr[i,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
        #=
        # Block high precision
        Q, R = bhh_QR(Ah, r);
        berr[i,t,1] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,1] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpBQR2
        Q, R = bhh_QR(A, r);
        berr[i,t,2] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,2] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpBQR3
        Q, R = mpbhh_QR(A, r);
        berr[i,t,3] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,3] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # TSQR high precision
        Q, R = par_TSQR(Ah, L);
        berr[i,t,4] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,4] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR2
        Q, R = par_TSQR(A, L);
        berr[i,t,5] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,5] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # mpTSQR3
        Q, R = mpTSQR(A, L);
        berr[i,t,6] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,6] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)

        # hhQr
        Q, R = hh_QR(Ah);
        berr[i,t,7] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad)
        oerr[i,t,7] = opnorm(Matrix{d}(Q')*Matrix{d}(Q)-I)
        =#
    end
    open("../txtfiles/"*name*"b.txt", "a") do io
        writedlm(io, [berr[:,t]'])
        #=writedlm(io, [berr[:,t,1]'])
        writedlm(io, [berr[:,t,2]'])
        writedlm(io, [berr[:,t,3]'])
        writedlm(io, [berr[:,t,4]'])
        writedlm(io, [berr[:,t,5]'])
        writedlm(io, [berr[:,t,6]'])
        writedlm(io, [berr[:,t,7]'])
        writedlm(io, ['\n'])
        =#
    end
    open("../txtfiles/"*name*"f.txt", "a") do io
        writedlm(io, [oerr[:,t]'])
        #=writedlm(io, [oerr[:,t,1]'])
        writedlm(io, [oerr[:,t,2]'])
        writedlm(io, [oerr[:,t,3]'])
        writedlm(io, [oerr[:,t,4]'])
        writedlm(io, [oerr[:,t,5]'])
        writedlm(io, [oerr[:,t,6]'])
        writedlm(io, [oerr[:,t,7]'])
        writedlm(io, ['\n'])
        =#
    end
end
