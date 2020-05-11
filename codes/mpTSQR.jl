using LinearAlgebra
#######################
## Define structure. ##
#######################

#=
Q's are stored as Q[1] and Q[2] where Q = vcat(Q[1], Q[2]). 
=#

# Householder vectors normalized to sqrt(2) s.t. 2.0/(v'v)=1.0 "n1"
struct par_TSQR_components_a
    A 
    V # V = Matrix containing householder vectors
    R 
    Q 
    par_TSQR_components_a(A::Array{T,2}, n_stop::Int=size(A)[2]) where T<:AbstractFloat = new(
        copy(A),                                         # A
        zeros(T, size(A)[1], n_stop),                    # V
        zeros(T, size(A)[2], size(A)[2]),                # R
        [zeros(T, floor(Int, size(A)[1]/2), size(A)[2]), # [Q[1], 
        zeros(T, ceil(Int, size(A)[1]/2), size(A)[2])],  # Q[2]]
        ) 
end

# Householder vectors normalized s.t. need to keep track of 2.0/(v'v)
struct par_TSQR_components_b
    A 
    b # Householder constants
    V # V = Matrix containing householder vectors
    R 
    Q 
    par_TSQR_components_b(A::Array{T,2}, n_stop::Int=size(A)[2]) where T<:AbstractFloat = new(
        copy(A),                                         # A
        zeros(T, n_stop),#size(A)[2]),                            # b
        zeros(T, size(A)[1], n_stop),                    # V
        zeros(T, size(A)[2], size(A)[2]),                # R
        [zeros(T, floor(Int, size(A)[1]/2), size(A)[2]), # [Q[1],
        zeros(T, ceil(Int, size(A)[1]/2), size(A)[2])],  # Q[2]]
        )
end

include("householderqr.jl")

#=
hhvec  = "n1": Each of the householder vectors are normalized to sqrt(2).
       = "n2": Each of the householder vectors are normalized such that the first entry is 1.
              If want_Q=false, return b vector that contains 2/(v'v). 
       = "nn": (Not normalized). Vectors calculated as is. 
              If want_Q=false, return b vector that contains 2/(v'v). 
=#

@inline function Vcolumns(h::Int, n::Int)
    if h == n 
        return  n-1
    else 
        return  n 
    end
end
######################################
## Initialize structure and problem ##
######################################
# no b vector #
function initProblem(A::Array{T,2}, L::Int; hhvec::String="n2") where T<:AbstractFloat
    m = size(A)[1]; n = size(A)[2];
    h = floor(Int, m/(2^L));           # height of small blocks
    ptc = Dict("n1"=> par_TSQR_components_a, "n2" => par_TSQR_components_b, "nn"=> par_TSQR_components_b)
    n_stop = Vcolumns(h, n)
    levels = Array{Array{ptc[hhvec],1},1}(undef, L+1)
    # Initialize number of blocks at each level.
    for i = 1 : length(levels)
        levels[i] = Array{ptc[hhvec],1}(undef, 2^(L-(i-1)))
    end
    # Initialize the first level by writing A as a block-column. 
    for j = 1 : 2^L
        if j == 2^L # Last block may have more rows. 
            n_stop = Vcolumns(m - (2^L-1)*h, n)
            levels[1][j] = ptc[hhvec](A[(j-1)*h+1: end, :], n_stop);
        else
            levels[1][j] = ptc[hhvec](A[(j-1)*h+1: j*h, :], n_stop);
        end
    end
    return levels
end

#############
## TSQR!!! ##
#############

@inline function binary_index(j::Int)
    k1 = ceil(Int, j/2)
    b = floor(Int, j/2)
    if k1 == b
        k2 = 2
    else
        k2 = 1
    end
    return k1, k2
end

function par_TSQR_compute_a(levels::Array{Array{par_TSQR_components_a,1},1}, A::Array{T,2}) where T<:AbstractFloat
    L = length(levels)-1
    m = size(A)[1]; n = size(A)[2];
    if L == 0
        return hh_QR(A, want_Q=true, pivot=false, thin=true, hhvec="n1")
    end
    # Forming R and storing householder vectors
    for i = 1 : L
        for j = 1 : 2^(L-i)
            levels[i][2*j-1].V[:, :], levels[i][2*j-1].R[:, :] = hh_QR(levels[i][2*j-1].A, want_Q=false, hhvec="n1")
            levels[i][2*j].V[:, :],   levels[i][2*j].R[:, :]   = hh_QR(levels[i][2*j].A,   want_Q=false, hhvec="n1")
            levels[i+1][j] = par_TSQR_components_a([levels[i][2*j-1].R; levels[i][2*j].R], n) # n_stop= n always since these A's are 2n by n. 
        end
    end
    levels[end][1] = par_TSQR_components_a([levels[L][1].R; levels[L][2].R], n);

    # Building Q.
    Q = Array{T,2}(undef, 0, n);
    tempQ, levels[end][1].R[:,:] = hh_QR(levels[end][1].A, want_Q=true, hhvec="n1");
    levels[end][1].Q[1][:,:] = tempQ[1:n, :]
    levels[end][1].Q[2][:,:] = tempQ[n+1:end, :]
    
    for i = L:-1:1
        for j = 1 : 2^(L-(i-1))
            k1, k2 = binary_index(j)
            tempQ = hh_multiply(
                levels[i][j].V, vcat(levels[i+1][k1].Q[k2], zeros(T, size(levels[i][j].A)[1]-size(levels[i+1][k1].Q[k2])[1], n)),
                [T(0)], 
                "n1"
                )
            if i == 1
                Q = vcat(Q, tempQ)
            else
                levels[i][j].Q[1] = tempQ[1: n, :]
                levels[i][j].Q[2] = tempQ[n+1: end, :]
            end
        end
    end
    return Q, UpperTriangular(levels[end][1].R)
end

function par_TSQR_compute_b(levels::Array{Array{par_TSQR_components_b,1},1}, A::Array{T,2}, hhvec::String="n2") where T<:AbstractFloat
    L = length(levels)-1
    m = size(A)[1]; n = size(A)[2];
    if L == 0
        return hh_QR(A, want_Q=true, pivot=false, thin=true, hhvec=hhvec)
    end
    # Forming R and storing householder vectors
    for i = 1 : L
        for j = 1 : 2^(L-i)
            levels[i][2*j-1].b[:], levels[i][2*j-1].V[:, :], levels[i][2*j-1].R[:, :] = hh_QR(levels[i][2*j-1].A, want_Q=false, hhvec=hhvec)
            levels[i][2*j].b[:],   levels[i][2*j].V[:, :],   levels[i][2*j].R[:, :]   = hh_QR(levels[i][2*j].A,   want_Q=false, hhvec=hhvec)
            levels[i+1][j] = par_TSQR_components_b([levels[i][2*j-1].R; levels[i][2*j].R], n) # n_stop= n always since these A's are 2n by n. 
        end
    end
    levels[end][1] = par_TSQR_components_b([levels[L][1].R; levels[L][2].R], n);

    # Building Q.
    Q = Array{T,2}(undef, 0, n);
    tempQ, levels[end][1].R[:,:] = hh_QR(levels[end][1].A, want_Q=true, hhvec=hhvec);
    levels[end][1].Q[1][:,:] = tempQ[1:n, :]
    levels[end][1].Q[2][:,:] = tempQ[n+1:end, :]
    
    for i = L:-1:1
        for j = 1 : 2^(L-(i-1))
            k1, k2 = binary_index(j)
            tempQ = hh_multiply(
                levels[i][j].V, 
                vcat(levels[i+1][k1].Q[k2], zeros(T, size(levels[i][j].A)[1]-size(levels[i+1][k1].Q[k2])[1], n)),
                levels[i][j].b, 
                hhvec
                )
            
            if i == 1
                Q = vcat(Q, tempQ)
            else
                levels[i][j].Q[1] = tempQ[1: n, :]
                levels[i][j].Q[2] = tempQ[n+1: end, :]
            end
        end
    end
    return Q, UpperTriangular(levels[end][1].R)
end


function par_TSQR(A::Array{T,2}, L::Int; hhvec::String="n2") where T<:AbstractFloat
    if hhvec ∉ ["n1", "n2", "nn"]
        error("Must choose n1, n2 or nn as hhqr normalization type.")
    else
        levels = initProblem(A, L, hhvec=hhvec)
    end
    if hhvec == "n1"
        Q, R = par_TSQR_compute_a(levels, A)
    else
        Q, R = par_TSQR_compute_b(levels,A, hhvec)
    end
    return Q, R
end

###############
# Diagnostics #
###############
using Printf

function numLevs!(A::Array{T,2}) where T<:AbstractFloat
    m,n=size(A);
    L = Int(floor(log2(m/n)));
    @printf("The maximum allowable number of levels is %d", L)
    return L
end