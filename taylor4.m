function y = taylor4(f, fp, fpp, fppp, h, n, t0, y0)
% TAYLOR4 Taylor series method assignment
%
% DESCRIPTION
%   y = TAYLOR4(f, fp, fpp, fppp, h, n, t0, y0) estimates y for t in
%   [t0, h*n] given it's first four derivates using Taylor's fourth order
%   method described as follows:
%       y(i+1) = y(i) + h*f(t(i),y(i)) + ... + h^4/4!*f'''(t(i),y(i))
%       t(i) = t0 + i*h
%       y(1) = y0
%
% INPUTS
%   f   - dy/dx(t,y)
%   fp  - d2y/dx2(t,y)
%   fpp - d3y/dx3(t,y)
%   fppp- d4y/dx4(t,y)
%   h   - step length in time
%   n   - number of subintervals in time
%   t0  - intitial time
%   y0  - y(t0)
%
% OUTPUTS
%   y   - vector of y for t in [t0, h*n]
%
% COPYRIGHT (C) Russell Maguire 2017

%% Time vector
t = t0+h*(0:n);

%% Preallocate solution
y = zeros(size(t));

% Apply intital condition
y(1) = y0;

%% Determine coefficients
c = h.^(0:4)./factorial(0:4);

%% Apply method
for i = 1:n
    y(i+1) = c * [y(i)
                  f(t(i),y(i))
                  fp(t(i),y(i))
                  fpp(t(i),y(i))
                  fppp(t(i),y(i))];
end
end