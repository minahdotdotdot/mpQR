using LinearAlgebra
l=Float16; h=Float32; d=Float64;

function genmat(m::Int,n::Int,κ::T,dt::DataType) where T<:AbstractFloat
	#α = (κ-1)/n;
	FL = qr(randn(h,m,n));
	FR = qr(randn(h,n,n));
	D = Diagonal(exp10.(-κ/(n-1)*(0:n-1)));
	#A = Matrix(F.Q)*(α*ones(d,n,n)+Matrix{d}(I, n,n))
	if dt == h
            return FL.Q*Matrix((FR.Q*D')')
        else
            return Matrix{dt}(FL.Q*Matrix((FR.Q*D')'))
        end
end
