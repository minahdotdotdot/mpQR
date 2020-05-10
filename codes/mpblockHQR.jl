include("householderqr.jl")
h = Float32
l = Float16

################
## WY update ### 
################

function buildWY(V::Matrix{h}, # Columns are householder vectors. (Should be lower trapezoidal.)
 b::Vector{h};                  # householder constants, if "n2" or "nn", just oens if "n1".
 )
	r = length(b); m̃ = size(V)[1];
	W = b[1]*V[1:end,1];
	for j = 2 : r
		z = b[j]*(V[:,j] - reshape(W,m̃,j-1) * (V[:,1:j-1]' *V[:,j]))
		W = hcat(W, z);      # is now m̃ by j
	end
	return W
end

function mpbhh_QR(
    A::Matrix{l}, r::Int;
    want_Q::Bool=true, 
    thin::Bool=true, 
    hhvec::String= "n2"
    )
	m,n = size(A);
	if m < n
	    error("System must be overdetermined or square.")
	end
	n_stop = n
	if m == n
	    n_stop = n-1
	end
    N = ceil(Int, n_stop/r);
    # Initialize Q, R, b, V, W, Y
    R  = deepcopy(A);
    λ = 1; k = 0;
    W = Vector{Any}(undef, N);
    Y = Vector{Any}(undef, N);
    while λ <= n_stop
        τ = min(λ+r-1, n); k += 1;
        b, Y[k], C =hh_QR_n2(Matrix{h}(R[λ:end, λ:τ]), want_Q=false, thin=false);
        W[k] = buildWY(Y[k], b);

        # Cast down.
        R[λ:end, λ:τ] = Matrix{l}(C);
        W[k] = Matrix{l}(W[k]);
        Y[k] = Matrix{l}(Y[k]);
        R[λ:end, τ+1:end] = mpWYupdate(W[k], Y[k], R[λ:end, τ+1:end])
        #R[λ:end, τ+1:end] -= Y[k]*(W[k]' * R[λ:end, τ+1:end]);
        λ = τ + 1;
    end
    if want_Q == true
	    if thin==true
	    	Q = Matrix{l}(I,m,n);
	    	R = R[1:n,1:n]
	       for i = N:-1:1
	    	  Q[(i-1)*r+1:end,(i-1)*r+1:end] -= W[i] * (Y[i]'*Q[(i-1)*r+1:end,(i-1)*r+1:end])
	       end
	       return Q, R
        end
    end
end

function mpWYupdate(W::Matrix{l}, Y::Matrix{l}, B::Matrix{l})
    h = Float32
    temp = bFMA(W, bFMA(-Y', B); C=B);
    return temp
end


# bFMA emulates TensorCore matrix-matrix multiply and accumulate
function bFMA(A::Matrix{l}, B::Matrix{l}; C=0)
    h=Float32;
    m,p=size(A); n=size(B)[2]
    if C == 0
        C = zeros(h,m,n);
    end
    return Matrix{l}(Matrix{h}(A)*Matrix{h}(B)+Matrix{h}(C))

end

function bFMA(A::Adjoint{l,Matrix{l}}, B::Matrix{l}; C=0)
    h=Float32;
    m,p=size(A); n=size(B)[2]
    if C == 0
        C = zeros(h,m,n);
    end
    return Matrix{l}(Matrix{h}(A)*Matrix{h}(B)+Matrix{h}(C))

end