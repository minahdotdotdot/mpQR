using DelimitedFiles, PyPlot,LaTeXStrings
l=Float16; h=Float32; d=Float64;
ul = .5*eps(l); uh = .5*eps(h); ud = .5*eps(d); 

function loadtxt(name::String, ave::Bool=true)
	b = readdlm("../txtfiles/"*name*"b.txt")[2:end,:];
	f = readdlm("../txtfiles/"*name*"f.txt")[2:end,:];
	if ave ==true
		return sum(b, dims=1)'/size(b)[1], sum(f, dims=1)'/size(f)[1]
	else
		return b, f
	end
end
cs = exp10.(0:0.25:2)   # A, D, G, J
cs2 = exp10.(0:0.5:4)   # A2, D2, G2, J2, Z
ns = 2 .^(7:11)         # B, H2
ns2 = 2 .^(7:12)        # B2
ns3 = 2 .^(6:10)        # H
ns4 = 2 .^(8:11)        # Y
rs = 2 .^(6:10)         # C, C2, I2, W
rs2 = 2 .^(5:9)         # I
ms = 2 .^(9:11)         # E, K
ms2 = 2 .^(10:13)       # E2, K2
ms3 = 2 .^(10:12)        # X
ms4 = 2 .^(8:13)        # X2
L = 1 : 5               # F, F2, L, L2

Yb, Yf = loadtxt("Y", false);
m=4096; r=128;
boundh = m*uh*ns4;
bound2 = ns4.^(1/2).*(ns4./r)*ul + ns4.^(3/2)*m/4*uh;
bound3 = ns4.^(3/2).*(100*ul + m/4*uh);
fig,ax = subplots()
#ax.set_yscale("log")
#ax.set_xscale("log")
#scatter(ns4, sum(Yb[7:8:end,:], dims=1)'/10, marker="x", s=200, label="hhQRh")
scatter(ns4, sum(Yb[1:8:end,:], dims=1)'/10, marker="x", s=200, label="BQRh", c=:blue)
#scatter(ns4, sum(Yb[4:8:end,:], dims=1)'/10, marker="x", s=200, label="TSQRh")

#scatter(ns4, sum(Yb[8:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="hhQRl")
#scatter(ns4, sum(Yb[2:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpBQR2", c=:red)
#scatter(ns4, sum(Yb[5:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpTSQR2")

scatter(ns4, sum(Yb[3:8:end,:], dims=1)'/10, alpha=0.5, label="mpBQR3", c=:green)
#scatter(ns4, sum(Yb[6:8:end,:], dims=1)'/10, alpha=0.5, label="mpTSQR3")
plot(ns4, boundh,label="BQR bound", c=:blue)
#plot(ns4, bound2,label="mpBQR2 bound", c=:red)
plot(ns4, bound3,label="mpBQR3 bound", c=:green)


axhline(.5*eps(l), label=L"u^{(l)}")
legend(bbox_to_anchor=(.9, .9))
xlim(minimum(ns4)*.9, maximum(ns4)/.9)
#ylim(minimum(Yb)*.9, maximum(Yb)/.9)
xlabel("Number of Columns")
ylabel("Frobenius norm Backward Error")
title("m=4096, c=1000, r=128")

#=
Wb, Wf = loadtxt("W", false);
fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
t = size(Wb)[1]/3;
m=4096;
n=4096;
scatter(rs, sum(Wb[1:3:end,:], dims=1)'/t, marker="x", s=200, label="BQRh")
scatter(rs, sum(Wb[2:3:end,:], dims=1)'/t, marker="D", alpha=0.5, label="mpBQR2")
#plot(rs, ((2 ./rs)*eps(Float16)/eps(Float32) .+m/4), label="mpBQR3 bound")
scatter(rs, sum(Wb[3:3:end,:], dims=1)'/t, alpha=0.5, label="mpBQR3")
legend(bbox_to_anchor=(.75, .5))
#xlim(minimum(rs)*.9, maximum(rs)/.9)
#ylim(minimum(Wb)*.5, 1)
xlabel("Size of Column block partition")
ylabel("Frobenius norm Backward Error")
title("Varying block size while m=n=4096")

fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(ms3, sum(Wb[7:8:end,:], dims=1)'/10, marker="x", s=200, label="hhQRh")
scatter(ms3, sum(Wb[1:8:end,:], dims=1)'/10, marker="x", s=200, label="BQRh")
scatter(ms3, sum(Wb[4:8:end,:], dims=1)'/10, marker="x", s=200, label="TSQRh")

scatter(ms3, sum(Wb[8:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="hhQRl")
scatter(ms3, sum(Wb[2:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpBQR2")
scatter(ms3, sum(Wb[5:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpTSQR2")

scatter(ms3, sum(Wb[3:8:end,:], dims=1)'/10, alpha=0.5, label="mpBQR3")
scatter(ms3, sum(Wb[6:8:end,:], dims=1)'/10, alpha=0.5, label="mpTSQR3")


legend(bbox_to_anchor=(.75, .5))
xlim(minimum(ms3)*.9, maximum(ms3)/.9)
xlabel("Number of Rows")
ylabel("Frobenius norm Backward Error")
title("n=256, r=256, L=2, c=1000")



Xb, Xf = loadtxt("X", false);
cs = exp10.(range(0, stop=4, length=9));
fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(ms3, sum(Xb[7:8:end,:], dims=1)'/10, marker="x", s=200, label="hhQRh")
scatter(ms3, sum(Xb[1:8:end,:], dims=1)'/10, marker="x", s=200, label="BQRh")
scatter(ms3, sum(Xb[4:8:end,:], dims=1)'/10, marker="x", s=200, label="TSQRh")

scatter(ms3, sum(Xb[8:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="hhQRl")
scatter(ms3, sum(Xb[2:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpBQR2")
scatter(ms3, sum(Xb[5:8:end,:], dims=1)'/10, marker="D", alpha=0.5, label="mpTSQR2")

scatter(ms3, sum(Xb[3:8:end,:], dims=1)'/10, alpha=0.5, label="mpBQR3")
scatter(ms3, sum(Xb[6:8:end,:], dims=1)'/10, alpha=0.5, label="mpTSQR3")


legend(bbox_to_anchor=(.75, .5))
xlim(minimum(ms3)*.9, maximum(ms3)/.9)
xlabel("Number of Rows")
ylabel("Frobenius norm Backward Error")
title("n=256, r=256, L=2, c=1000")
=#


#=
Ab, Af = loadtxt("A2")
mpAb, mpAf = loadtxt("mpA2")
#rs = 2 .^(6:10);
cs = exp10.(range(0, stop=4, length=9));
#semilogy(cs, Ab, label="A-b")
#semilogy(cs, Af, label="A-f")
loglog(cs, mpAb, label="mpA-b")
#semilogy(cs, mpAf, label="mpA-f")
xlabel("Condition Numbers")
ylabel("2-norm Relative error")
title("BQR3")
legend()
=#


#=
Ab, Af = loadtxt("A")
cs = exp10.(range(0, stop=2, length=9));
semilogy(cs, sum(Ab, dims=1)'/size(Ab)[1], label="backward")
semilogy(cs, sum(Af, dims=1)'/size(Af)[1], label="orthog")
xlabel("Condition Numbers")
ylabel("2-norm Relative error")
legend()

Bb, Bf = loadtxt("B")
ns = 2 .^(7:11)
loglog(ns, sum(Bb, dims=1)'/size(Bb)[1], label="backward")
loglog(ns, sum(Bf, dims=1)'/size(Bf)[1], label="orthog")
xlabel("# of Columns")
ylabel("2-norm Relative error")
legend()


Cb, Cf = loadtxt("C")
rs = 2 .^(6:10);
loglog(rs, sum(Cb, dims=1)'/size(Cb)[1], label="backward")
loglog(rs, sum(Cf, dims=1)'/size(Cf)[1], label="orthog")
xlabel("Column Partition size")
ylabel("2-norm Relative error")
legend()
=#