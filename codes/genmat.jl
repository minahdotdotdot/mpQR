using LinearAlgebra
l=Float16; h=Float32; d=Float64;

function genmat(m::Int,n::Int,κ::T,dt::DataType) where T<:AbstractFloat
	α = (κ-1)/n;
	F = qr(randn(m,n));
	A = Matrix(F.Q)*(α*ones(d,n,n)+Matrix{d}(I, n,n))
	return Matrix{dt}(A/norm(A))
end
