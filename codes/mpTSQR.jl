using LinearAlgebra
l=Float16; h=Float32; d=Float64;

#######################
## Define structure. ##
#######################

#=
Q's are stored as Q[1] and Q[2] where Q = vcat(Q[1], Q[2]). 
=#
struct par_TSQR_components
    A 
    b
    V
    W
    R
    Q
    par_TSQR_components(A::Matrix{h}, n_stop::Int=size(A)[2]) = new(
        copy(A),                                      # A
        zeros(h, n_stop),                                # b  
        zeros(l, size(A)[1], n_stop),                    # V
        zeros(l, size(A)[1], n_stop),                    # W
        zeros(h, size(A)[2], size(A)[2]),                # R
        [zeros(l, floor(Int, size(A)[1]/2), size(A)[2]), # [Q[1],
        zeros(l, ceil(Int, size(A)[1]/2), size(A)[2])],  # Q[2]]
        )
end

include("householderqr.jl")
include("mp.jl")

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
function initProblem(A::Matrix{h}, L::Int)
    m = size(A)[1]; n = size(A)[2];
    hh = floor(Int, m/(2^L));           # height of small blocks
    n_stop = Vcolumns(hh, n)
    levels = Array{Array{par_TSQR_components,1},1}(undef, L+1)
    # Initialize number of blocks at each level.
    for i = 1 : length(levels)
        levels[i] = Array{par_TSQR_components,1}(undef, 2^(L-(i-1)))
    end
    # Initialize the first level by writing A as a block-column and stored in h. 
    for j = 1 : 2^L
        if j == 2^L # Last block may have more rows. 
            n_stop = Vcolumns(m - (2^L-1)*hh, n)
            levels[1][j] = par_TSQR_components(Matrix{h}(A[(j-1)*hh+1: end, :]), n_stop);
        else
            levels[1][j] = par_TSQR_components(Matrix{h}(A[(j-1)*hh+1: j*hh, :]), n_stop);
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

function par_TSQR_compute(levels::Array{Array{par_TSQR_components,1},1}, A::Matrix{h}, hhvec::String="n2")
    L = length(levels)-1
    m = size(A)[1]; n = size(A)[2];
    if L == 0
        return hh_QR(A, want_Q=true, pivot=false, thin=true, hhvec=hhvec)
    end
    # Forming R and storing householder vectors
    for i = 1 : L
        for j = 1 : 2^(L-i)
            levels[i][2*j-1].b[:], V1, levels[i][2*j-1].R[:,:] = hh_QR(levels[i][2*j-1].A, want_Q=false, hhvec=hhvec)
            levels[i][2*j].b[:], V2,   levels[i][2*j].R[:,:]   = hh_QR(levels[i][2*j].A,   want_Q=false, hhvec=hhvec)
            levels[i+1][j] = par_TSQR_components([levels[i][2*j-1].R; levels[i][2*j].R], n) # n_stop= n always since these A's are 2n by n. 

            #Build WY.
            W1 = buildWY(V1, levels[i][2*j-1].b[:])
            W2 = buildWY(V2, levels[i][2*j].b[:])

            #Cast down  (they can only be stored in l)
            levels[i][2*j-1].V[:,:] = Matrix{l}(V1)
            levels[i][2*j].V[:,:]   = Matrix{l}(V2)
            levels[i][2*j-1].W[:,:] = Matrix{l}(W1)
            levels[i][2*j].W[:,:]   = Matrix{l}(W2)
        end
    end
    #end = L+1
    levels[end][1] = par_TSQR_components([levels[L][1].R; levels[L][2].R], n);

    # Building Q.
    Q = Matrix{l}(undef, 0, n);
    tempQ, levels[end][1].R[:,:] = hh_QR(levels[end][1].A, want_Q=true, hhvec=hhvec);

    #Cast down
    levels[end][1].R[:,:] = Matrix{l}(levels[end][1].R[:,:])
    tempQ = Matrix{l}(tempQ);

    levels[end][1].Q[1][:,:] = tempQ[1:n, :]
    levels[end][1].Q[2][:,:] = tempQ[n+1:end, :]
    
    for i = L:-1:1
        for j = 1 : 2^(L-(i-1))

            #Figure out sizes
            k1, k2 = binary_index(j)
            m̃ = size(levels[i][j].A)[1]-size(levels[i+1][k1].Q[k2])[1];

            #Apply Q to a factor from level below.
            tempQ = mpWYupdate(levels[i][j].W, levels[i][j].V, vcat(levels[i+1][k1].Q[k2], zeros(l, m̃,n)))
            #bFMA(levels[i][j].Q, vcat(levels[i+1][k1].QQ[k2], zeros(l, m̃,n) ))
            if i == 1
                # If at the topmost level, concatenate the 2^L blocks to form Q. 
                Q = vcat(Q, tempQ)
            else
                #At each level, bisect the Q factor so it can be multiplied in the level above.
                levels[i][j].Q[1] = tempQ[1: n, :]
                levels[i][j].Q[2] = tempQ[n+1: end, :]
            end
        end
    end
    return Q, UpperTriangular(levels[end][1].R)
end

function mpTSQR(A::Matrix{l}, L::Int; hhvec::String="n2")
    Ah = Matrix{h}(A);
    if hhvec ∉ ["n1", "n2", "nn"]
        error("Must choose n1, n2 or nn as hhqr normalization type.")
    else
        levels = initProblem(Ah, L)
    end
    Q, R = par_TSQR_compute(levels, Ah, hhvec)
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