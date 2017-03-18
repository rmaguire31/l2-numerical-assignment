function [x,t,u] = conv_diff_eig(a,b,n,T,m,c,v,d,f)
% CONV_DIFF Convection-diffusion-reaction FDM assignment with eigenvalues
%
% DESCRIPTION
%   Computes the particular solution for a convection-diffusion-reaction
%   problem in 1D with a proportional reaction component:
%       u_t(x,t) = c^2*u_xx(x,t) + v*u_x(x,t) + d*u(x,t) in [a,b]x[0,T]
%       BC: u(a,t) = u(b,t) = 0
%       IC: u(x,0) = f(x)
%
%   This problem is approximated using the following formulae:
%       u_t(x,t) ~ DF_k u(x,t)|x = (u(x,t+k) - u(x,t))/k
%       u_x(x,t) ~ DF_h u(x,t)|t = (u(x+h,t) - u(x,t))/h
%       u_xx(x,t) ~ D2C_h u(x,t)|t = (u(x+h,t) - 2*u(x,t) + u(x-h,t))/h^2
%
%   Which reduces to the problem to the following stencil:
%                  [ d*k + 1 ] [0  1  0][u(i-1,j)]
%       u(i,j+1) = [  v*k/h  ].[0 -1  1][ u(i,j) ]
%                  [c^2*k/h^2] [1 -2  1][u(i+1,1)]
%
%   The stencil is used to produce a tridiagonal matrix for finding the
%   solution u(x,j+1) from u(x,j) at time j, starting from u(x,0) which is
%   given.
%                             [s1 s2  0 ...  0]
%                             [s0 s1 s2 ...  0]
%       u(x,j+1) = S*u(x,j) = [ 0 s0 s1 ...  0] u(x,j)
%                             [ :  :  :  '.  :]
%                             [ 0  0  0 ... s1]
%
%             [s0]   [0  1  0][ d*k + 1 ]
%       where [s1] = [0 -1  1][  v*k/h  ]
%             [s2]   [1 -2  1][c^2*k/h^2]

%   By diagonalising and elementwise-exponentiating a block diagonal matrix
%   is constructed for finding the solution u(x,t) for all t in (0,T] from
%   u(x,0) which is given.
%                [S   0 ...   0][f(x)]
%       u(x,t) = [0 S^2 ...   0][f(x)] for t in (0,t]
%                [:   :  '.   0][  : ]
%                [0   0 ... S^m][f(x)]
% INPUTS
%   a -  lower limit for interval in x
%   b -  upper limit for interval in x
%   n -  number of subintervals in x
%   T -  upper limit for interval in t
%   m -  number of steps in t
%   c -  diffusion coefficient
%   v -  convection coefficient
%   d -  reaction coefficient
%   f -  function describing x at initial time t = 0
%
% OUTPUTS
%   x -  space vector in [a,b]
%   t -  time vector in [0,T]
%   u -  grid of solutions in [a,b]x[0,T]
%
% COPYRIGHT (C) Russell Maguire 2017

% Step lengths
h = (b-a)/n;
k = T/m;

% Time and space vector
x = a:h:b;
t = 0:k:T;

% Construct stencil
s = [d*k+1 v*k/h c^2*k/h^2] * [0  1 0   % Zeroth order
                               0 -1 1   % First order
                               1 -2 1]; % Second order
S = cell(m, 1);
S{1} = spdiags(repmat(s, n-1, 1), -1:1, n-1, n-1);

[V,D1,W] = eig(full(S{1}));
D = D1;
for j = 1:m-1
    D = D .* D1;
    S{j} = sparse(real(V*D*W));
end
S = blkdiag(S{:});

% Preallocate solution
u = zeros(length(x),length(t));

% Apply intial condition
u(2:n,1) = f(x(2:n));

% Apply Boundary conditions
u([1 end], :) = 0;

% Apply stencil
v = real(S * repmat(u(2:n,1), m, 1));
u(2:n,2:end) = reshape(v, n-1, []);
end