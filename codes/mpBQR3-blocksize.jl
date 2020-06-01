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

#rs3 = 2 .^(2:12)        # dat
#m=2^13; n=m;

name = "1106_"
rs=ceil.(Int,2 .^range(1, stop=8, length=19))
rlen = length(rs);
trials=2;
m=2^11; n=2^6;
include("W.jl")
datb, datf = loadtxt(name, false);

datbh = reshape(datb[:,1],rlen,trials); datbh = maximum(datbh, dims=2);
datbmp3 = reshape(datb[:,2],rlen,trials); datbmp3 = maximum(datbmp3, dims=2);
datfh = reshape(datf[:,1],rlen,trials); datfh = maximum(datfh, dims=2);
datfmp3 = reshape(datf[:,2],rlen,trials); datfmp3 = maximum(datfmp3, dims=2);
#datbh = datb[1:2:end,1];
#datbmp3 = datb[2:2:end,1];
#datfh = datf[1:2:end,1];
#datfmp3 = datf[2:2:end,1];

boundh = n^(3/2)*m*uh*ones(19)#length(datbh));
bound3 = n^(3/2)*(10 ./rs*ul .+ m/4*uh);

fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(rs, datbh, label="BQRh backward", c=:blue)
scatter(rs, datbmp3, label="mpBQR3 backward", c=:green)
scatter(rs, datfh, label="BQRh orthogonal", c=:blue, marker="x")
scatter(rs, datfmp3, label="mpBQR3 orthogonal", c=:green, marker="x")

plot(rs, boundh,label=L"n^{3/2}mu^{(h)}", c=:blue)
plot(rs, bound3,label=L"n^{1/2}(10Nu^{(l)}+\frac{nm}{4}u^{(h)})", c=:green)
#plot(rs3, boundh,label="BQR bound", c=:blue)
#plot(rs3, bound2,label="mpBQR2 bound", c=:red)
#plot(rs3, bound3,label="mpBQR3 bound", c=:green)
axhline(.5*eps(h), label=L"u^{(fp32)}", c=:black)
axhline(.5*eps(l), label=L"u^{(fp16)}", c=:black, linestyle=:dashed)

xticks()
xlim(minimum(rs)*.9, maximum(rs)/.9)
ylim(.1*eps(h), maximum(bound3)/.9)
xlabel("r: block size")
ylabel("Norm error")
title("BQRh and mpBQR3 performance on "*string(m)*"-by-"*string(n)*" matrices")
legend(bbox_to_anchor=(0.1,0.22,.75,.2), ncol=2)
#legend(bbox_to_anchor=(.15, .5), ncol=2)