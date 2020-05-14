using DelimitedFiles, PyPlot

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
rs = 2 .^(6:10)         # C, C2, I2
rs2 = 2 .^(5:9)         # I, 
ms = 2 .^(9:11)         # E, K
ms2 = 2 .^(10:13)       # E2, K2
L = 1 : 5               # F, F2, L, L2

Xb, Xf = loadtxt("X", false);
cs = exp10.(range(0, stop=4, length=9));
fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(ms, sum(Xb[1:8:end,:], dims=1)'/5, marker="o", s=200, label="BQRh")
#scatter(ms, sum(Xb[2:8:end,:], dims=1)'/5, marker="D", alpha=0.5, label="mpBQR2")
scatter(ms, sum(Xb[3:8:end,:], dims=1)'/5, alpha=0.5, label="mpBQR3")
scatter(ms, sum(Xb[4:8:end,:], dims=1)'/5, marker="x", s=200, label="TSQRh")
#scatter(ms, sum(Xb[5:8:end,:], dims=1)'/5, marker="D", alpha=0.5, label="mpTSQR2")
scatter(ms, sum(Xb[6:8:end,:], dims=1)'/5, alpha=0.5, label="mpTSQR3")
scatter(ms, sum(Xb[7:8:end,:], dims=1)'/5, marker="x", s=200, label="hhQRh")

scatter(ms, sum(Xb[1:8:end,:], dims=1)'/5, marker="x", s=200, label="BQRh")
#scatter(ms, sum(Xb[7:8:end,:], dims=1)'/5, marker="D", alpha=0.5, label="hhQRl")
legend(bbox_to_anchor=(.75, .5))
xlim(minimum(ms)*.9, maximum(ms)/.9)
xlabel("Number of Rows")
ylabel("Frobenius norm Backward Error")
title("n=256, r=256, L=2, c=1000")


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