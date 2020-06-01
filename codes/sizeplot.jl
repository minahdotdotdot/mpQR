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

ms = exp10.(range(3, stop=log10(45000), length=14));
ns = floor.(Int,ms/4);
ms = ceil.(Int, ms);
L = 2;
κs = 1e3*ones(length(ms));
dt = h;
sb, sf = loadtxt("size_", false);
ssb=Vector{d}(sb[1:9*floor(Int,size(sb)[1]/9),1]);
ssf=Vector{d}(sf[1:9*floor(Int,size(sf)[1]/9),1]);
ssb = reshape(ssb,(9,Int(length(ssb)/9)))
ssf = reshape(ssf,(9,Int(length(ssf)/9)))
#=
boundh = n^(3/2)*m/4*uh*ones(length(W2bh));
bound3 = n^(3/2)*(25 ./rs3*ul .+ m/4*uh);
=#

fig,ax = subplots()
ax.set_yscale("log")
ax.set_xscale("log")
scatter(ms[1:size(ssb)[2]], ssb[1:9:end,:]', marker="+", c=:red  , label="mpHQR2") 
scatter(ms[1:size(ssb)[2]], ssb[3:9:end,:]', marker="+", c=:blue , label="mpBQR2") 
scatter(ms[1:size(ssb)[2]], ssb[6:9:end,:]', marker="+", c=:green, label="mpTSQR2") 

scatter(ms[1:size(ssb)[2]], ssb[4:9:end,:]', marker="D", c=:blue, label="mpBQR3") 
scatter(ms[1:size(ssb)[2]], ssb[7:9:end,:]', marker="D", c=:green, label="mpTSQR3") 

scatter(ms[1:size(ssb)[2]], ssb[8:9:end,:]', marker="x", c=:red, label="HQRh") 
scatter(ms[1:size(ssb)[2]], ssb[2:9:end,:]', marker="x", c=:blue, label="BQRh") 
scatter(ms[1:size(ssb)[2]], ssb[5:9:end,:]', marker="x", c=:green, label="TSQRh") 

axhline(.5*eps(h), label=L"u^{(fp32)}", c=:black)
axhline(.5*eps(l), label=L"u^{(fp16)}", c=:black, linestyle=:dashed)
axhline(1, label="1",c=:red, linewidth=0.5)
#=
plot(rs3, boundh,label="BQR bound", c=:blue)
#plot(rs3, bound2,label="mpBQR2 bound", c=:red)
plot(rs3, bound3,label="mpBQR3 bound", c=:green)
=#
xticks()
xlim(minimum(ms)*.9, maximum(ms)/.9)
ylim(.1*eps(h), maximum(ssb)/.8)
xlabel("m, number of rows")
ylabel("Backward Frobenius norm error")
title(L"QR relative backward errors of $m$-by-$\frac{m}{4}$ matrices")
legend()
#legend(bbox_to_anchor=(.95, .8))
