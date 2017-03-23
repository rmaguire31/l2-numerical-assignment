function y = taylor4(f, fp, fpp, fppp, h, n, t0, y0)
% TAYLOR4 Taylor series method assignment
%
% DESCRIPTION
%   y = TAYLOR4(f, fp, fpp, fppp, h, n, t0, y0) estimates y for t in
%   [t0, h*n] given it's first four derivates using Taylor's fourth order
%   method described as follows:
%       y_i+1 = y_i + h*f(t_i,y_i) + ... + h^4/4!*f'''(t_i,y_i)
%       t_i+1 = t_i + h
%       given y_0, t_0 and h
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

%% Determine coefficients
c = h.^(0:4)./factorial(0:4);

%% Apply initial conditions
y = y0;
t = t0;

%% Apply method
for i = 1:n
    y = c * [y
             f(t,y)
             fp(t,y)
             fpp(t,y)
             fppp(t,y)];
    t = t + h;
end
end