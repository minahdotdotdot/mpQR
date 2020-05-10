using LinearAlgebra, Random
#= 
First option is default. 

want_Q = true: Returns an explicit Q matrix 
       = false: Returns matrix V whose columns are householder vectors

pivot  = false: only returns Q, R. 
       = true: Returns permutation vector P, Q, R, where P_1 = sparse(P, collect(1:n), ones(n)). AP_1 = QR
       

thin   = true: If want_Q=true, gives, H_1...H_nI_{m\times n}
       = false: If want_Q=true, gives, H_1...H_nI_{m\times m}

hhvec  = "n2": Each of the householder vectors are normalized such that the first entry is 1.
              If want_Q=false, return b vector that contains 2/(v'v). 
       = "n1": Each of the householder vectors are normalized to sqrt(2).
       = "nn": (Not normalized). Vectors calculated as is. 
              If want_Q=false, return b vector that contains 2/(v'v). 
=#

###########
## hhvec ## 
###########
@inline function hh_vec(x::Array{T,1}) where T<:AbstractFloat
    num_type = eltype(x)
    v = x
    ss = 1
    if x[1] != 0
        ss = sign(x[1])
    end
    sigma = -ss * norm(x,2)                # sigma = -sign(x[1])*\|x\|_2
    v[1] = x[1]-sigma                      # v = x - sigma*e_1
    v = num_type(1.0)/sqrt(-sigma*v[1])*v  # v normalized s.t. v'v = 2.
    return sigma, v
end

@inline function hh_vec_n2(x::Array{T,1}) where T<:AbstractFloat # Default option: "n2"
    num_type = eltype(x)
    v = x
    ss = 1
    if x[1] != 0
        ss = sign(x[1])
    end
    sigma = -ss * norm(x,2)                # sigma = -sign(x[1])*\|x\|_2
    v[1] = x[1]-sigma                      # v = x - sigma*e_1
    bbeta = -v[1]/sigma; # calculation of bbeta simplified by algebra. (less dot products!)
    return sigma, bbeta, v/v[1]
end

@inline function hh_vec_nn(x::Array{T,1}) where T<:AbstractFloat
    num_type = eltype(x)
    v = x
    ss = 1
    if x[1] != 0
        ss = sign(x[1])
    end
    sigma = -ss * norm(x,2)                # sigma = -sign(x[1])*\|x\|_2
    v[1] = x[1]-sigma                      # v = x - sigma*e_1  
    #bbeta = num_type(1.0)/(-sigma*v[1]) # calculation of bbeta simplified by algebra. (less dot products!)
    return sigma, num_type(1.0)/(-sigma*v[1]), v
end

################
## hhmultiply ## 
################
@inline function hh_multiply(
	V::Array{T,2},        # Columns are householder vectors. (Should be lower trapezoidal.)
	Q::Array{T,2},        # Matrix we are multiplying the householder reflections to. 
	b::Array{T,1},        # householder constants, if "n2" or "nn". 
	hhvec::String = "n2"  # type or normalization on householder vectors (or no normalization.)
	) where T<:AbstractFloat
	if size(V)[1] != size(Q)[1]
		error("V and Q must have the same number of rows.")
	end
	n_stop = size(V)[2]
	if size(V)[1] == size(V)[2]
		n_stop -=1
	end
	if hhvec == "n1"
		for i = n_stop : -1: 1
        	Q[i:end, i:end] = Q[i:end, i:end] - V[i:end, i]*(V[i:end, i]'*Q[i:end, i:end])
        end
    elseif hhvec =="n2" || hhvec=="nn"
	    for i = n_stop : -1: 1
	        Q[i:end, i:end] = Q[i:end, i:end] - b[i]*V[i:end, i]*(V[i:end, i]'*Q[i:end, i:end])
	    end
    else
        error("hhvec needs to be n1, n2, or nn.\n")
	end
    return Q 
end

###########
## hhQR ## 
###########

function hh_QR_nn(
    A::Array{T,2}; 
    want_Q::Bool=true, 
    pivot::Bool=false, 
    thin::Bool=true
    ) where T<:AbstractFloat
    m = size(A)[1]; n = size(A)[2]
    if m < n
        error("System must be overdetermined or square.")
    end
    n_stop = n
    if m == n
        n_stop = n-1
    end

    # Initialize  b, V, R, P.
    b = ones(T, n_stop) # if hhvec=n1, this will be eliminated in hh_multiply.
    V = zeros(T, m, n_stop)
    R  = deepcopy(A)
    if pivot == true
        P = collect(1:n)
    end

    for i = 1 : n_stop
        # Column pivot.
        if pivot == true
            # Permute [i:end] columns, and leave [1:i-1] alone. 
            tempP = [collect(1:i-1); sortperm(vec(sum(R[i:m, i:n].^2, dims=1)), rev=true) .+(i-1)]
            P = P[tempP]    # Update permutation. 
            R = R[:, tempP] # Update R via permutation.
        end

        # Do only if entries below diagonal are not all zero.        
        if maximum(abs.(R[i+1:m, i])) > T(0)
            # Save householder vector, v and constant, beta. 
            sigma, b[i], V[i:m, i]= hh_vec_nn(R[i:m, i])

            # Update R. 
            R[i, i]    = sigma                                                            # newR[i,i] = -sigma
            R[i+1:m, i] = zeros(T, m-i, 1)                                                # newR ith column below ith row is zeros.
            R[i:m, i+1:n] = R[i:m, i+1:n] - b[i]*V[i:m, i]*(V[i:m, i]'*R[i:m, i+1:n]) # newR, i:m rows by i+1:n columns are updated.
        end
    end

    if thin ==true
            Q = Matrix{T}(I, m, n)
            R = UpperTriangular(R[1:n, :])
    else
            Q = Matrix{T}(I, m, m)
    end

    if want_Q == false
        if pivot == true
            return b, V, R, P
        else
            return b, V, R
        end
    else
        #Q = hh_multiply(V, Q, b, hhvec)
        if pivot == true
            return hh_multiply(V, Q, b, "nn"), R, P
        else 
            return hh_multiply(V, Q, b, "nn"), R
        end
    end
end


function hh_QR_n2(
    A::Array{T,2}; 
    want_Q::Bool=true, 
    pivot::Bool=false, 
    thin::Bool=true
    ) where T<:AbstractFloat
    m = size(A)[1]; n = size(A)[2]
    if m < n
        error("System must be overdetermined or square.")
    end
    n_stop = n
    if m == n
        n_stop = n-1
    end

    # Initialize  b, V, R, P.
    b = ones(T, n_stop) # if hhvec=n1, this will be eliminated in hh_multiply.
    V = zeros(T, m, n_stop)
    R  = deepcopy(A)
    if pivot == true
        P = collect(1:n)
    end

    for i = 1 : n_stop
        # Column pivot.
        if pivot == true
            # Permute [i:end] columns, and leave [1:i-1] alone. 
            tempP = [collect(1:i-1); sortperm(vec(sum(R[i:m, i:n].^2, dims=1)), rev=true) .+ (i-1)]
            P = P[tempP]    # Update permutation. 
            R = R[:, tempP] # Update R via permutation.
        end

        # Do only if entries below diagonal are not all zero.        
        if maximum(abs.(R[i+1:m, i])) > T(0)
            # Save householder vector, v and constant, beta. 
            sigma, b[i], V[i:m, i]= hh_vec_n2(R[i:m, i])    


            # Update R. 
            R[i, i]    = sigma                                                        # newR[i,i] = -sigma
            R[i+1:m, i] = zeros(T, m-i, 1)                                            # newR ith column below ith row is zeros.
            R[i:m, i+1:n] = R[i:m, i+1:n] - b[i]*V[i:m, i]*(V[i:m, i]'*R[i:m, i+1:n]) # newR, i:m rows by i+1:n columns are updated. 
        end
    end

    if thin ==true
            Q = Matrix{T}(I, m, n)
            R = UpperTriangular(R[1:n, :])
    else
            Q = Matrix{T}(I, m, m)
    end

    if want_Q == false
        if pivot == true
            return b, V, R, P
        else
            return b, V, R
        end
    else
        #Q = hh_multiply(V, Q, b, hhvec)
        if pivot == true
            return hh_multiply(V, Q, b, "n2"), R, P
        else 
            return hh_multiply(V, Q, b, "n2"), R
        end
    end
end

function hh_QR(
    A::Array{T,2}; 
    want_Q::Bool=true, 
    pivot::Bool=false, 
    thin::Bool=true, 
    hhvec::String= "n2"
    ) where T<:AbstractFloat
    if hhvec == "n2"
        return hh_QR_n2(A, want_Q=want_Q, pivot=pivot, thin=thin)
    elseif hhvec =="nn"
        return hh_QR_nn(A, want_Q=want_Q, pivot=pivot, thin=thin)
    elseif hhvec == "n1"
        m = size(A)[1]; n = size(A)[2]
        if m < n
            error("System must be overdetermined or square.")
        end
        n_stop = n
        if m == n
            n_stop = n-1
        end

        # Initialize  V, R, P.
        V = zeros(T, m, n_stop)
        R  = deepcopy(A)
        if pivot == true
            P = collect(1:n)
        end

        for i = 1 : n_stop
            # Column pivot.
            if pivot == true
                # Permute [i:end] columns, and leave [1:i-1] alone. 
                tempP = [collect(1:i-1); sortperm(vec(sum(R[i:m, i:n].^2, dims1)), rev=true) .+(i-1)]
                P = P[tempP]    # Update permutation. 
                R = R[:, tempP] # Update R via permutation.
            end

            # Do only if entries below diagonal are not all zero.        
            if maximum(abs.(R[i+1:m, i])) > T(0)
                # Save householder vector, v and constant, beta. 
                sigma, V[i:m, i] = hh_vec(R[i:m, i])

                # Update R. 
                R[i, i]    = sigma                                                            # newR[i,i] = -sigma
                R[i+1:m, i] = zeros(T, m-i, 1)                                                # newR ith column below ith row is zeros.
                R[i:m, i+1:n] = R[i:m, i+1:n] -V[i:m, i]*(V[i:m, i]'*R[i:m, i+1:n])      # newR, i:m rows by i+1:n columns are updated.
            end
        end

        if thin ==true
                Q = Matrix{T}(I, m, n)
                R = UpperTriangular(R[1:n, :])
        else
                Q = Matrix{T}(I, m, m)
        end

        if want_Q == false
            if pivot == true
                return V, R, P
            else
                return V, R
            end
        else
            #Q = hh_multiply(V, Q, b, hhvec)
            if pivot == true
                return hh_multiply(V, Q, [T(1)], hhvec), R, P
            else
                return hh_multiply(V, Q, [T(1)], hhvec), R
            end
        end
    else
        error("hhvec needs to be n1, n2, or nn.\n")
    end
end


@inline function dottest(x::Array{T,1}, y::Array{T,1}) where T<: Float16
    S = Array{Float32,1}
    z = x .* y
    return abs(sum(S(z)) -  S(x)'*S(y))
end

function backsolve(R::UpperTriangular{T,Array{T,2}}, b::Array{T,1}) where T<:AbstractFloat
    PR = eltype(R)
    if size(R)[2] != length(b)
        error("Sizes for R and b do not match!")
    else
        x = zeros(PR, length(b))
        for i = length(b) : -1 : 1
            x[i] = (b[i] - R[i, i+1:end]'*x[i+1:end]) / R[i,i]
        end
        return x
    end
end

#=
trials = 10000;
T=Float16; 
n = [100, 1000, 10000, 100000, 1000000]
stats = zeros(Float32,length(n),2)
for j = 1 : length(n)
    x = zeros(Float32, trials)
    for i = 1 : trials
        x[i] = dottest(randn(T,n[j]), randn(T,n[j]))
    end
    stats[j,:] = [mean(x),var(x)]
end

using PyPlot
fig, ax = subplots()
ax[:set_yscale]("log")
ax[:set_xscale]("log")
scatter(n, stats[:,1], label="mean")
scatter(n, stats[:,2], label="var")
legend()
=#

#=
@inline function hh_vec{T<:AbstractFloat}(x::Array{T,1}, hhvec::String = "n1")
    num_type = eltype(x)
    v = x
    ss = 1
    if x[1] != 0
        ss = sign(x[1])
    end
    sigma = -ss * norm(x,2)                # sigma = -sign(x[1])*\|x\|_2
    v[1] = x[1]-sigma                      # v = x - sigma*e_1
    
    if hhvec == "n2" 
        bbeta = num_type(1.0)/(-sigma*v[1]); # calculation of bbeta simplified by algebra. (less dot products!)
        return sigma, bbeta*v[1]^2, v/v[1]
    elseif hhvec == "nn"
        #bbeta = num_type(1.0)/(-sigma*v[1]) # calculation of bbeta simplified by algebra. (less dot products!)
        return sigma, num_type(1.0)/(-sigma*v[1]), v
    elseif hhvec == "n1"
        v = num_type(1.0)/sqrt(-sigma*v[1])*v  # v normalized s.t. v'v = 2.
        return sigma, v
    else
        error("hhvec needs to be n1, n2, or nn.\n")
    end
end

function hh_QR{T<:AbstractFloat}(
    A::Array{T,2}; 
    want_Q::Bool=true, 
    pivot::Bool=false, 
    thin::Bool=true, 
    hhvec::String= "n1"
    )
    m = size(A)[1]; n = size(A)[2]
    if m < n
        error("System must be overdetermined or square.")
    end
    if hhvec ∉ ["n1", "n2", "nn"]
        error("Must choose n1, n2 or nn as hhqr normalization type.")
    end
    n_stop = n
    if m == n
        n_stop = n-1
    end

    # Initialize  b, V, R, P.
    b = ones(T, n_stop) # if hhvec=n1, this will be eliminated in hh_multiply.
    V = zeros(T, m, n_stop)
    R  = deepcopy(A)
    if pivot == true
        P = collect(1:n)
    end

    for i = 1 : n_stop
        # Column pivot.
        if pivot == true
            # Permute [i:end] columns, and leave [1:i-1] alone. 
            tempP = [collect(1:i-1); i - 1 + sortperm(vec(sum(R[i:m, i:n].^2, 1)), rev=true)]
            P = P[tempP]    # Update permutation. 
            R = R[:, tempP] # Update R via permutation.
        end

        # Do only if entries below diagonal are not all zero.        
        if maximum(abs.(R[i+1:m, i])) > T(0)
            # Save householder vector, v and constant, beta. 
            if hhvec == "n2" || hhvec =="nn"
                sigma, b[i], V[i:m, i]= hh_vec(R[i:m, i], hhvec)    
            else
                sigma, V[i:m, i]= hh_vec(R[i:m, i], hhvec)
            end

            # Update R. 
            R[i, i]    = sigma                                                            # newR[i,i] = -sigma
            R[i+1:m, i] = zeros(T, m-i, 1)                                                # newR ith column below ith row is zeros.
            if hhvec == "n2" || hhvec =="nn"
                R[i:m, i+1:n] = R[i:m, i+1:n] - b[i]*V[i:m, i]*(V[i:m, i]'*R[i:m, i+1:n]) # newR, i:m rows by i+1:n columns are updated.     
            else
                R[i:m, i+1:n] = R[i:m, i+1:n] -V[i:m, i]*(V[i:m, i]'*R[i:m, i+1:n])      # newR, i:m rows by i+1:n columns are updated.  
            end
        end
    end

    if thin ==true
            Q = Matrix{T}(I, m, n)
            R = R[1:n, :]
    else
            Q = Matrix{T}(I, m, m)
    end

    if want_Q == false
        if pivot == true
            if hhvec == "n1"
                return V, R, P
            else
                return b, V, R, P
            end
        else
            if hhvec == "n1"
                return V, R
            else
                return b, V, R
            end
        end
    else
        #Q = hh_multiply(V, Q, b, hhvec)
        if pivot == true
            if hhvec == "n1"
                return hh_multiply(V, Q, [T(1)], hhvec), R, P
            else
                return hh_multiply(V, Q, b, hhvec), R, P
            end
        else 
            if hhvec == "n1"
                return hh_multiply(V, Q, [T(1)], hhvec), R
            else
                return hh_multiply(V, Q, b, hhvec), R
            end
        end
    end
end



=#
