function gamma(u,m)
	return (m*u)/(1-(m*u))
end

ms = ceil.(Int,range(1_000,stop=700_000,length=7_000));
ns = ceil.(Int,range(100,stop=40_000,length=1_000));

HQRerr=zeros()
TSQRL2err=
TSQRL4er=
for (i,m) in enumerate(ms)
	for (j,n) in enumerate(ns)
		
	end
end