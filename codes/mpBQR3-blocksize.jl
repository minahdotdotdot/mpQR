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

rs3 = 2 .^(2:12)        # W2
m=2^13; n=m;
W2b, W2f = loadtxt("W2", false);
W2bh = W2b[1:2:end,1];
W2bmp3 = W2b[2:2:end,1];
W2fh = W2f[1:2:end,1];
W2fmp3 = W2f[2:2:end,1];

boundh = n^(3/2)*m/4*uh*ones(length(W2bh));
bound3 = n^(3/2)*(25 ./rs3*ul .+ m/4*uh);

fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(rs3, W2bh, label="BQRh backward", c=:blue)
scatter(rs3, W2bmp3, label="mpBQR3 backward", c=:green)
scatter(rs3, W2fh, label="BQRh orthogonal", c=:blue, marker="x")
scatter(rs3, W2fmp3, label="mpBQR3 orthogonal", c=:green, marker="x")

plot(rs3, boundh,label="BQR bound", c=:blue)
#plot(rs3, bound2,label="mpBQR2 bound", c=:red)
plot(rs3, bound3,label="mpBQR3 bound", c=:green)
axhline(.5*eps(h), label=L"u^{(fp32)}", c=:black)
axhline(.5*eps(l), label=L"u^{(fp16)}", c=:black, linestyle=:dashed)

xticks()
xlim(minimum(rs3)*.9, maximum(rs3)/.9)
ylim(.1*eps(h), maximum(bound3)/.9)
xlabel("r: block size")
ylabel("Norm error")
title("BQRh and mpBQR3 performance on 8192-by-8192 matrices")
legend(bbox_to_anchor=(.15, .5), ncol=2)