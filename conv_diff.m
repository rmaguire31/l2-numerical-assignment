function [x,t,u] = conv_diff(a,b,n,T,m,c,v,d,f)
% CONV_DIFF Convection-diffusion-reaction FDM assignment
%
% DESCRIPTION
%   [x,t,u] = CONV_DIFF(a,b,n,T,m,c,v,d,f) computes the particular
%   solution for a convection-diffusion-reaction problem in 1D with a
%   linear reaction component:
%       u_t(x,t) = c^2*u_xx(x,t) + v*u_x(x,t) + d*u(x,t) in [a,b]x[0,T]
%       BC: u(a,t) = u(b,t) = 0
%       IC: u(x,0) = f(x)
%
%   This problem is approximated using the following formulae:
%       u_t(x,t) ~ DF_k u(x,t)|x = (u(x,t+k) - u(x,t))/k
%       u_x(x,t) ~ DF_h u(x,t)|t = (u(x+h,t) - u(x,t))/h
%       u_xx(x,t) ~ D2C_h u(x,t)|t = (u(x+h,t) - 2*u(x,t) + u(x-h,t))/h^2
%   NB: Forward difference is chosen over central difference for it's
%   increased stability.
%
%   Which reduces to the problem to the following stencil:
%                  [ d*k + 1 ] [ 0  1  0][u(i-1,j)]
%       u(i,j+1) = [0.5*v*k/h].[ 0 -1  1][ u(i,j) ]
%                  [c^2*k/h^2] [ 1 -2  1][u(i+1,j)]
%
%   The stencil is used to produce a tridiagonal matrix for finding the
%   solution u(x,j+1) from u(x,j) at time j, starting from u(x,0) which is
%   given.
%                  [s1 s2  0 ...  0  0]
%                  [s0 s1 s2 ...  0  0]
%       u(x,j+1) = [ 0 s0 s1 ...  0  0] u(x,j)
%                  [ :  :  :  '.  :  :]
%                  [ 0  0  0 ... s1 s2]
%                  [ 0  0  0 ... s2 s1]
%
%             [s0]   [ 0  1  0][ d*k + 1 ]
%       where [s1] = [ 0 -1  1][  v*k/h  ]
%             [s2]   [ 1 -2  1][c^2*k/h^2]
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

%% Space and time vectors
h = (b-a)/n;
k = T/m;

x = a:h:b;
t = 0:k:T;

%% Preallocate solution
u = zeros(length(x),length(t));

% Apply intial condition
u(2:n,1) = f(x(2:n));

% Apply Boundary conditions
u([1 end], :) = 0;

%% Construct stencil for each time slot.
% The zero boundary conditions mean rows and columns for u_0 and u_n can be
% removed giving an order n-1 stencil.
s = [d*k+1 v*k/h c^2*k/h^2]*[ 0  1  0
                              0 -1  1    % DC
                              1 -2  1 ]; % D2C
S = spdiags(repmat(s, n-1, 1), -1:1, n-1, n-1);

%% Apply stencil at each time 
for j = 1:m
    u(2:n,j+1) = S * u(2:n,j);
end
end