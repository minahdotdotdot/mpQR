
using LinearAlgebra
function gammaf(m, u)
	return m*u/(1-(m*u))
end

ms = range(1_000, stop=100_000, length=10_000);
ns = range(100, stop=5_000, length=500);

errs=zeros(3,10_000,500)

for (i,m) in enumerate(ms)
	for (j,n) in enumerate(ns)
		if m >= n
			errs[1,i,j] = n^(3/2)*gammaf(m,.5*eps(Float64))
			if m >= 4*n
				errs[2,i,j] = n^(3/2)*(gammaf(m*2^(-2),.5*eps(Float64))+2*gammaf(2*n,.5*eps(Float64)))
				if m >= 16*n
					errs[3,i,j] = n^(3/2)*(gammaf(m*2^(-4),.5*eps(Float64))+4*gammaf(2*n,.5*eps(Float64)))
				end
			end
		end
	end
end

titles=["", "L=2", "L=4"];
using PyPlot
PyPlot.matplotlib.rc("font",size=12)
Nr = 1
Nc = 3
images = [];
fig, axs = plt.subplots(Nr, Nc,figsize=(6,5))

fig.suptitle("HQR                TSQR with L levels      ")
vmin=minimum(log10.(errs[errs.!=0]))
vmax=maximum(log10.(errs[:]))
#vmin=minimum(errs[errs.!=0])
#vmax=maximum(errs[:])

for i = 1 : Nc
	#push!(images,axs[i].imshow(errs[i,:,:],
	push!(images,axs[i].imshow(log10.(errs[i,:,:]),
		cmap="gist_rainbow",
		extent=(ns[1],ns[end],ms[1],ms[end]),
		vmin=vmin,vmax=vmax,origin="lower left",
		aspect=2*(ns[end]-ns[1])/(ms[end]-ms[1])#,
		#interpolation="spline16"
		))
	if i==1
		axs[i].set_ylabel("Number of rows (m)")
		axs[i].ticklabel_format(scilimits=(0,4))
		#axs[i].set_yticks(zSPREAD[:])
		#axs[i].set_yticklabels(string.(zSPREAD)[:])
	else
		#axs[i].set_yticks([])
		axs[i].set_yticklabels([])
	end
	if i ==2
		axs[i].set_xlabel("Number of columns (n)")
	end
	axs[i].set_title(titles[i])
end
#=fig.subplots_adjust(top=0.88,
bottom=0.11,
left=0.11,
right=0.9,
hspace=0.2,
wspace=0.2)=#
#fig.tight_layout()
fig.colorbar(images[1], ax=axs, orientation="horizontal", fraction=.05)
