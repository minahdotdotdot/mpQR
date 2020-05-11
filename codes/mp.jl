function mpWYupdate(W::Matrix{l}, Y::Matrix{l}, B::Matrix{l};b=4)
    h = Float32
    temp = bFMA(W, bFMA(-Y', B,b=b); C=B,b=b);
    return temp
end


# bFMA emulates TensorCore matrix-matrix multiply and accumulate in mixed precision
function bFMA(A::T, B::Matrix{l}; C=0, b=4) where T<: Union{Matrix{l}, Adjoint{l,Matrix{l}}}
    h=Float32;
    #return Matrix{l}(Matrix{h}(A)*Matrix{h}(B)+Matrix{h}(C))
    return Matrix{l}(bMMM(Matrix{h}(A), Matrix{h}(B),C=C, b=b))
end


# bMMM emulates block matrix-matrix multiply and accumulate in uniform precision
function bMMM(A::Matrix{T}, B::Matrix{T}; C=0, b::Int) where T<:AbstractFloat
    m,p=size(A); n=size(B)[2];
    mm = ceil(Int, m/b)*b; pp = ceil(Int, p/b)*b; nn = ceil(Int, n/b)*b;
    AA = zeros(T, mm, pp); AA[1:m,1:p]=A;
    BB = zeros(T, pp, nn); BB[1:p,1:n]=B;
    if C == 0
        CC = zeros(T,mm,nn);
    else
        CC = zeros(T, mm, nn); CC[1:m,1:n]=Matrix{T}(C);
    end
    x=0
    while x < mm
        y=0
        while y < nn
            z=0
            while z < pp
                CC[x+1:x+b, y+1:y+b] += AA[x+1:x+b, z+1:z+b]*BB[z+1:z+b, y+1:y+b];
                z+=b
            end
            y+=b
        end
        x+=b
    end
    return CC[1:m,1:n]
end


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