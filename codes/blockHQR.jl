include("householderqr.jl")

################
## WY update ### 
################

function WYupdate(V::Matrix{T}, # Columns are householder vectors. (Should be lower trapezoidal.)
 b::Vector{T};                  # householder constants, if "n2" or "nn", just oens if "n1".
 ) where T<:AbstractFloat
	r = length(b); m̃ = size(V)[1];
	W = b[1]*V[1:end,1];
	for j = 2 : r
		z = b[j]*(V[:,j] - reshape(W,m̃,j-1) * (V[:,1:j-1]' *V[:,j]))
		W = hcat(W, z);      # is now m̃ by j
	end
	return W
end

function bhh_QR_n2(
    A::Array{T,2}, r::Int, m::Int, n::Int, n_stop::Int, N::Int;
    want_Q::Bool=true, 
    thin::Bool=true
    ) where T<:AbstractFloat
    # Initialize Q, R, b, V, W, Y
    R  = deepcopy(A);
    λ = 1; k = 0;
    W = Vector{Matrix{T}}(undef, N);
    Y = Vector{Matrix{T}}(undef, N);
    while λ <= n_stop
    	τ = min(λ+r-1, n); k += 1;
    	b, Y[k], R[λ:end, λ:τ] =hh_QR_n2(R[λ:end, λ:τ], want_Q=false, thin=false);
    	W[k] = WYupdate(Y[k], b);
    	R[λ:end, τ+1:end] -= Y[k]*(W[k]' * R[λ:end, τ+1:end]);
    	λ = τ + 1;
    end
    return W,Y,R
end

function bhh_QR_nn(
    A::Array{T,2}, r::Int, m::Int, n::Int, n_stop::Int, N::Int;
    want_Q::Bool=true, 
    thin::Bool=true
    ) where T<:AbstractFloat
    # Initialize Q, R, b, V, W, Y
    R  = deepcopy(A);
    λ = 1; k = 0;
    N = ceil(Int, n_stop/r);
    W = Vector{Matrix{T}}(undef, N);
    Y = Vector{Matrix{T}}(undef, N);
    while λ <= n_stop
    	τ = min(λ+r-1, n); k += 1;
    	b, Y[k], R[λ:end, λ:τ] =hh_QR_nn(R[λ:end, λ:τ], want_Q=false, thin=false);
    	W[k]= WYupdate(Y[k], b);
    	R[λ:end, τ+1:end] -= Y[k]*(W[k]' * R[λ:end, τ+1:end]);
    	λ = τ + 1;
    end
    return W,Y,R
end

function bhh_QR(
    A::Matrix{T}, r::Int;
    want_Q::Bool=true, 
    thin::Bool=true, 
    hhvec::String= "n2"
    ) where T<:AbstractFloat
	m,n = size(A);
	if m < n
	    error("System must be overdetermined or square.")
	end
	n_stop = n
	if m == n
	    n_stop = n-1
	end
    N = ceil(Int, n_stop/r);
    if hhvec == "n2"
        W,Y,R = bhh_QR_n2(A, r, m, n, n_stop, N, want_Q=want_Q, thin=thin)
    elseif hhvec =="nn"
        W,Y,R= bhh_QR_nn(A, r, m, n, n_stop, N, want_Q=want_Q, thin=thin)
    elseif hhvec == "n1"
	    # Initialize Q, R, b, V, W, Y
	    R  = deepcopy(A);
	    λ = 1; k = 0;
	    W = Vector{Matrix{T}}(undef, N);
	    Y = Vector{Matrix{T}}(undef, N);
	    while λ <= n_stop
	    	τ = min(λ+r-1, n); k += 1;
	    	Y[k], R[λ:end, λ:τ] =hh_QR_n1(R[λ:end, λ:τ], want_Q=false, thin=false);
	    	W[k] = WYupdate(Y[k], ones(T,τ-λ+1));
	    	R[λ:end, τ+1:end] -= Y[k]*(W[k]' * R[λ:end, τ+1:end]);
	    	λ = τ + 1;
	    end
	end
    if want_Q == true
	    if thin==true
	    	Q = Matrix{T}(I,m,n);
	    	R = R[1:n,1:n]
	    else
	    	Q = Matrix{T}(I,m,m)
	    end
	    for i = N:-1:1
	    	Q[(i-1)*r+1:end,(i-1)*r+1:end] -= W[i] * (Y[i]'*Q[(i-1)*r+1:end,(i-1)*r+1:end])
	    end
	    return Q, R
	else
		if thin == true
			return W, Y, R[1:n, 1:n]
		else
			return W,Y,R
		end
	end
end