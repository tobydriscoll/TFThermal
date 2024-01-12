module Chebyshev

using LinearAlgebra, FFTW

export points, chebeval, diffmat,  interpmat, interp, clencurt

function points(n,interval=[-1.,1.])
    #return [ -cos( pi*k/(n-1) ) for k = 0:n-1 ];
    x = [ sin(pi*k/(2(n-1))) for k = 1-n:2:n-1 ]
    return (x.+1)/2*(interval[2]-interval[1]) .+ interval[1]
end

"""
    chebpoly(n[,x])

Evaluate degree-`n` Chebyshev polynomial at given `x`. If `x` is not given, return a
callable function for the evaluation.
"""
function chebpoly(n::Integer,x::Number)
    return cos( n*acos(x) )
end

function chebpoly(n::Integer)
    return x -> chebpoly(n,x)
end

"""
    chebeval(c,x)

Evaluate at `x` the polynomial whose Chebyshev coefficients are given in Vector `c`
in increasing order. Uses Clenshaw's algorithm.
"""
function chebeval(c,x::Float64)
    bk2 = bk1 = 0.0;
    for k = length(c):-1:2
        #bk = c[k] + 2x*bk1 - bk2;
        bk = muladd(2x,bk1,c[k]-bk2);
        bk2 = bk1;
        bk1 = bk;
    end
    return muladd(x,bk1,c[1]-bk2);
end

function chebevalgrid(C,x,y)
    s = size(C);
    V = zeros(length(x),s[2]);
    for j = 1:s[2]
        c = view(C,:,j);
        for i = eachindex(x)
            V[i,j] = chebeval(c,x[i]);
        end
    end
    U = zeros(length(x),length(y));
    for i = 1:length(x)
        v = V[i,:];
        for j = eachindex(y)
            U[i,j] = chebeval(v,y[j]);
        end
    end
    return U
end

function chebevalgrid(C,x...)
    assert(ndims(C)==length(x));
    for d = 1:ndims(C)
        s = setindex!([size(C)...],length(x[d]),d)
        U = zeros(s...)
        for jpre in CartesianRange((s[1:d-1]...,))
            for jpost in CartesianRange((s[d+1:end]...,))
                j = (jpre,:,jpost);
                c = view(C,j...);
                for i = eachindex(x[d])
                    U[jpre,i,jpost] = chebeval(c,x[d][i]);
                end
            end
        end
        C = U;
    end
    return C
end

# Adapted from Differentiation Matrix Suite by Weideman and Reddy

function chebdif(N, M)
     ondiag = 1:N+1:N^2                   # Linear indices of diagonal elements

     n1,n2 = Int(floor(N/2)), Int(ceil(N/2))        # Indices used for flipping trick.
     th = [ k*pi/(2(N-1)) for k = 0:N-1 ]
     x = points(N)

     DX = 2*sin.(th'.+th).*sin.(th'.-th)          # Trigonometric identity.
     DX = [ DX[1:n1,:]; -reverse(reverse(DX[1:n2,:],dims=1),dims=2) ];   # Flipping trick.
     DX[ondiag] .= 1.0

     C = [ (-1.0)^(i+j) for i = 0:N-1, j=0:N-1 ]
     C[[1,N],:] *= 2
     C[:,[1,N]] /= 2

     Z = 1 ./ DX                         # Z contains entries 1/(x(k)-x(j))
     Z[ondiag] .= 0                      # with zeros on the diagonal.

     D = diagm(0=>ones(N))               # D contains diff. matrices.
     DM = ()
     for ell = 1:M
         D = ell*Z.*(C.*repmat(diag(D),1,N) - D); # Off-diagonals
         D[ondiag] = -sum(D,dims=2)                            # Correct main diagonal
         DM = (DM...,D);
     end

     return DM
end

chebdif(N) = chebdif(N,1)[1]

function diffmat(n,interval=[-1.,1.])
  x = -points(n);
  #n = n-1;
  n2 = ceil(Int,n/2);
  #D  = zeros(n+1,n2);
  D = zeros(n,n)

  # First column is special
  D[1,1] = -(2*(n-1)^2 + 1)/6;
  r = 2:n-1;
  @. D[r,1] = (-1)^r / (2*(x[r]-1));
  D[n,1] = (-1)^(n+1) / 2;


  # remainder of first n/2 cols
  for j = 2:n2
      for i = 1:n
          if i==j
              D[i,j] = x[j]/(2*(1-x[j]^2));
          else
              D[i,j] = (-1)^(i+j+1) / (x[i]-x[j]);
          end
          if (i==1 || i==n)
              D[i,j] *= 2
          end
      end
  end

  # Skew-flip symmetry
  for j = 1:n2
      for i = 1:n
          if iseven(n) || j<n2
            D[n+1-i,n+1-j] = -D[i,j]
            end
      end
  end

  # negative sum trick
  D[1:n+1:end] .= 0;
  D[1:n+1:end] = -sum(D,dims=2);
  return D*(2/(interval[2]-interval[1]))
end

"""
    interp(fvals,x)

Evaluates at `x` the interpolant whose values at 2nd-kind points are `fvals`.
"""
function interp(y,x::Number)
    n = length(y);
    xk = points(n);

    if n == 1
      return y[1]
    end

    numer = denom = 0.0;
    for k = 1:n
        a = (-1)^(n+k)/(x-xk[k]);
        if isinf(a)
            return y[k]
        end
        k==1 || k==n ? a/=2 : nothing;
        numer = muladd(y[k],a,numer);
        denom += a;
    end
    fx = numer/denom;
    return fx
end

function interp(Y,x...)
    assert(ndims(Y)==length(x));
    for d = 1:ndims(Y)
        s = setindex!([size(Y)...],length(x[d]),d);
        U = zeros(s...)
        for jpre in CartesianRange((s[1:d-1]...,))
            for jpost in CartesianRange((s[d+1:end]...,))
                j = (jpre,:,jpost)
                y = view(Y,j...)
                for i = eachindex(x[d])
                    U[jpre,i,jpost] = interp(y,x[d][i]);
                end
            end
        end
        Y = U;
    end
    return Y
end

function interpmat(y,x)

  N = length(x);
  M = length(y);

  # Nothing to do here!
  if ( M == N && all(x == y) )
      return diagm(0=>ones(N))
  end

  # Default to the Chebyshev barycentric weights:
  w = ones(N);
  w[2:2:end] .= -1;
  w[[1, N]] = 0.5*w[[1, N]];
  # Ensure w is a row vector:
  w = w';

  # Repmat(Y-X'):
  B = y .- x';

  # Construct the matrix:
  B = w./B;  # implicitly repeated down rows
  c = 1 ./ sum(B, dims=2);
  B = B.*c;

  # Where points coincide there will be division by zeros (as with bary.m).
  # Replace these entries with the identity:
  B[isnan.(B)] .= 1;

  return B
end

function coeffs(y)
    n = size(y)
    c = FFTW.r2r(y,FFTW.REDFT00)/prod(n.-1)
    # Adjust for the "forward x" direction
    for i in Iterators.product(2:2:k for k in n)
        c[i...] .= -c[i...]
    end
    return c
end

function clencurt(N)
    N = N-1;
    theta = pi*(0:N)/N;
    x = cos.(theta);
    w = zeros(N+1); ii = 2:N;
    v = ones(N-1);
    if mod(N,2)==0
        w[1] = 1/(N^2-1); w[N+1] = w[1];
        for k=1:N/2-1
            v -= 2*cos.(2*k*theta[ii])/(4*k^2-1);
        end
        v -= cos.(N*theta[ii])/(N^2-1);
    else
        w[1] = 1/N^2; w[N+1] = w[1];
        for k=1:(N-1)/2
            v -= 2*cos.(2*k*theta[ii])/(4*k^2-1);
        end
    end
    w[ii] = 2*v/N;
    return w
end


end # module
