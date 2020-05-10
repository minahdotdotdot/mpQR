using PyPlot, Printf
include("blockHQR.jl")
include("mpblockHQR.jl")

l = Float16; h=Float32; d=Float64;
m = 5000;
ns = [100, 250, 500, 1000, 2500, 5000]
trials = 10;
berr = zeros(d, length(ns), 3, trials);
ferr = zeros(d, length(ns), 3, trials);
r = 250;
for (i,n) in enumerate(ns)
        @printf("\n n=%d, trial: ",n)
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
	    berr[i,2,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)
	    ferr[i,2,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)


	    # High precision 
	    Q, R = bhh_QR(Ah, r);
	    berr[i,3,t] = norm(Matrix{d}(Q)*Matrix{d}(R)-Ad,2)
	    ferr[i,3,t] = norm(Matrix{d}(Q')*Matrix{d}(Q)-I,2)
        end
end



#=
m=100; p = 50; n = 30
A = randn(l, m, p); Ah = Matrix{h}(A);
B = randn(l, p, n); Bh = Matrix{h}(B);
C = randn(l, m, n); Ch = Matrix{h}(C);
F = randn(l, n, n); Fh = Matrix{h}(F)
D = bFMA(bFMA(A,B,C=C),F)
Dh = (Ah*Bh+Ch)*Fh;
Dnaive = (A*B+C)*Fh;
#Dip = bFMA(A, B, C=C, ip=true);
vmin = minimum(vcat(D,Dh,Dnaive)-vcat(Dh,Dnaive,D))
vmax = maximum(vcat(D,Dh,Dnaive)-vcat(Dh,Dnaive,D))

subplot(131)
imshow(Dh-D,vmin=vmin, vmax=vmax, cmap="seismic")
title("f32 vs bFMA")
#imshow(Dh-Dip,vmin=vmin, vmax=vmax, cmap="seismic")
#title("f32 vs bFMAip")
subplot(132)
imshow(Dh-Dnaive,vmin=vmin, vmax=vmax, cmap="seismic")
title("f32 vs julia f16")
subplot(133)
colorbar(imshow(Matrix{h}(D-Dnaive),vmin=vmin, vmax=vmax, cmap="seismic"))
title("bFMA vs julia f16")
#colorbar(imshow(Matrix{h}(Dip-Dnaive),vmin=vmin, vmax=vmax, cmap="seismic"))
#title("bFMAip vs julia f16")
=#
