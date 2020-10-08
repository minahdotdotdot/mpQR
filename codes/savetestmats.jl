include("genmat.jl")
using MAT
function savetestmats!(ms::Vector{Int}, ns::Vector{Int}, κs,
	dt::DataType, name::String)
	for i = 1 : length(ms)
		varname=name*string(i,pad=2)
		file = matopen("../"*varname*".mat", "w")
		write(file, varname, genmat(ms[i],ns[i],κs[i],dt))
		close(file)
	end
end
#start=11;
ms = exp10.(range(3, stop=log10(14000), length=10));
ns = 1000. *ones(10);
ms = ceil.(Int, ms);
κs = 1e3*ones(length(ms));
dt = h;
name="st3"
savetestmats!(ms, ns, κs, dt, name)
#savetestmats!(ms[1:1], ns[1:1], κs[1:1], dt, name)




