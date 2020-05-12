using DelimitedFiles, PyPlot

function loadtxt(name::String)
	b = readdlm("../txtfiles/"*name*"b.txt");
	f = readdlm("../txtfiles/"*name*"f.txt");
	return b[2:end,:],f[2:end,:]
end

#=
Ab, Af = loadtxt("A")
cs = exp10.(range(0, stop=2, length=9));
semilogy(cs, sum(Ab, dims=1)'/size(Ab)[1], label="backward")
semilogy(cs, sum(Af, dims=1)'/size(Af)[1], label="orthog")
xlabel("Condition Numbers")
ylabel("2-norm Relative error")
legend()
=#

#=
Bb, Bf = loadtxt("B")
ns = 2 .^(7:11)
loglog(ns, sum(Bb, dims=1)'/size(Bb)[1], label="backward")
loglog(ns, sum(Bf, dims=1)'/size(Bf)[1], label="orthog")
xlabel("# of Columns")
ylabel("2-norm Relative error")
legend()
=#

Cb, Cf = loadtxt("C")
rs = 2 .^(6:10);
loglog(rs, sum(Cb, dims=1)'/size(Cb)[1], label="backward")
loglog(rs, sum(Cf, dims=1)'/size(Cf)[1], label="orthog")
xlabel("Column Partition size")
ylabel("2-norm Relative error")
legend()