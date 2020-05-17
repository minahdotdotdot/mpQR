include("genmat.jl")
using MAT
ms = exp10.(range(3, stop=log10(45000), length=14));
ns = floor.(Int,ms/4);
ms = ceil.(Int, ms);
κs = 1e3*ones(length(ms));
dt = h;
name="W3"



function savetestmats!(ms::Vector{Int}, ns::Vector{Int}, κs,
	dt::DataType, name::String)
	for i = 1 : length(ms)
		varname=name*string(i,pad=2)
		file = matopen("../"*varname*".mat", "w")
		write(file, varname, genmat(ms[i],ns[i],κs[i],dt))
		close(file)
	end
end

savetestmats!(ms, ns, κs, dt, name)




