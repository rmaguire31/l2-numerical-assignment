function [x,t,u] = conv_diff_jordan(a,b,n,T,m,c,v,d,f)
% CONV_DIFF Convection-diffusion-reaction FDM assignment with jordan-normal
%
% DESCRIPTION
%   Computes the particular solution for a convection-diffusion-reaction
%   problem in 1D with a linear reaction component:
%       u_t(x,t) = c^2*u_xx(x,t) + v*u_x(x,t) + d*u(x,t) in [a,b]x[0,T]
%       u(a,t) = u(b,t) = 0 (Boundary Conditions)
%       u(x,0) = f(x) (Initial Condition)
%   This problem is approximated using the following formulae:
%       u_t(x,t) = (u(x,t+k) - u(x,t))/k (Forward Difference)
%       u_x(x,t) = (u(x+h,t) - u(x,t))/h (Forward Difference)
%       u_xx(x,t) = (u(x+h,t) - 2*u(x,t) + u(x-h,t))/h^2
%                                           (2nd Order Central Difference)
%   Which reduces to the following stencil:
%       u(i,j+1) = s0*u(i-1,j) + s1*u(i,j) + s2*u(i+1,j)
%       s0 = c^2*k/h^2
%       s1 = -2*c^2*k/h^2 - v*k/h + d*k + 1
%       s2 = c^2*k/h^2 + v*k/h
%   The stencil is used to produce a tridiagonal matrix for finding the
%   solution u(x,j+1) from u(x,j) at time j, starting from u(x,0) which is
%   given.
%       By diagonalising and elementwise-exponentiating a block diagonal
%   matrix is constructed for finding the solution u(x,t) for all t in
%   (0,T] from u(x,0) which is given. The Jordan-normal form of the
%   tridiagonal matrix is used rather than the eigenvalues, so some
%   numerical error is introduced.
%
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

[V,J1,W] = jordan(full(S{1}));
J = J1;
for j = 1:m-1
    J = J .* J1;
    S{j} = sparse(real(V*J*W));
end
S = blkdiag(S{:});

% Preallocate solution
u = zeros(length(x),length(t));

% Apply intial condition
u(2:n,1) = f(x(2:n));

% Apply Boundary conditions
u([1 end], :) = 0;

% Apply stencil
% for j = 1:m
%     u(2:n,j+1) = S * u(2:n,j);
% end
v = real(S * repmat(u(2:n,1), m, 1));
u(2:n,2:end) = reshape(v, n-1, []);
end