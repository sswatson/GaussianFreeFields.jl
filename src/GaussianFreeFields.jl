module GaussianFreeFields

export DGFF,
       fix_boundary_values

function DGFF(n::Int)
    h = complex(zeros(n,n));
    for j=1:n
        for k=1:n
            h[j,k] = 1/sqrt(2.0)*(j+k == 2 ? 0 : (randn() + im*randn()) * 1/sqrt(sin((j-1)*pi/n)^2+ sin((k-1)*pi/n)^2));
        end
    end
    return real(n*1/sqrt(2)*ifft(h))
end

function fix_boundary_values(h::Array{Float64,2},boundary_values::Array{Float64,2}=zeros(size(h)[1],size(h)[2]))
    n = size(h)[1]
    A = sparse([1,n],[1,n],[1,1.0])
    B = speye(n)
    C = B - A
    D_dense = zeros(n,n); 
    for i=1:n-1 
        D_dense[i,i+1]=-1; 
        D_dense[i+1,i]=-1; 
        D_dense[i,i] = 4; 
    end
    D_dense[1,1] = 1
    D_dense[1,2] = 0 
    D_dense[n,n-1] = 0
    D_dense[n,n] = 1
    D = sparse(D_dense);
    E_dense = zeros(n,n)
    for i=1:n-1 
        E_dense[i,i+1]=-1
        E_dense[i+1,i]=-1 
    end
    E_dense[n,n-1] = 0
    E_dense[1,2] = 0
    E = sparse(E_dense)
    M = kron(A,B) + kron(C,D) + kron(E,C)
    boundary = zeros(n^2)
    for i in 1:n-1
        boundary[i] = h[1,i] - boundary_values[1,i]
        boundary[n^2-n+i] = h[n,i] - boundary_values[n,i]
        boundary[n*i+1] = h[i+1,1] - boundary_values[i+1,1]
        boundary[n*i] = h[i,n] - boundary_values[i,n]
    end
    boundary[n^2] = h[n,n]
    
    return h - transpose(reshape(full(lufact(M) \ boundary),n,n))
end


end # module
