using LinearAlgebra, DelimitedFiles, MAT
l=Float16;d=Float64;h=Float32
function halfdot(x::Vector{l}, y::Vector{l})
    s = x[1]*y[1]
    for i = 2 : length(x)
        s+=x[i]*y[i]
    end
    return s
end
  
trials = 2*10^6;
m=1024;
errs = Vector{d}(undef,0);
for i = 1:trials
    x = rand(l,m);xd = Vector{d}(x);
    y = rand(l,m);yd = Vector{d}(y);
    err = abs(dot(xd,yd)-halfdot(x,y))/dot(abs.(xd),abs.(yd))
    #print(err,"\n")
    push!(errs, err)
end

function stats(errs::Vector{d})
    m=mean(errs)
    return [m, stdm(errs,m),maximum(errs)];
end

ustats = stats(errs);

errs = Vector{d}(undef,0);
for i = 1:trials
    x = randn(l,m);xd = Vector{d}(x);
    y = randn(l,m);yd = Vector{d}(y);
    err = abs(dot(xd,yd)-halfdot(x,y))/dot(abs.(xd),abs.(yd))
    #print(err,"\n")
    push!(errs, err)
end
nstats = stats(errs);


    
