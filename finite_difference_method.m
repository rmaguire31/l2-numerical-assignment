function finite_difference_method()
% FINITE_DIFFERENCE_METHOD Calculate the values for Table 3
%
% COPYRIGHT (C) Russell Maguire 2017

xc = 0.5;
tc = 0.2;

a = 0;
b = 1;
T = 0.2;
c = 1;
v = 0.5;
d = 0.02;
f = @(x) 4*x - 4*x.^2;

n = [4 16 32];
m = [10 100 200 400 800];

% Cartesian product of m and n.
[M, N] = meshgrid(m, n);
mn = [M(:) N(:)];

table3 = cell(size(mn,1), 4);
for i = 1:size(mn)
    m = mn(i,1);
    n = mn(i,2);
    [x, t, u] = conv_diff(a, b, n, T, m, c, v, d, f);
    
    table3(i,1) = {sprintf('%.3g', n)};
    table3(i,2) = {sprintf('%.3g', m)};
    table3(i,3) = {sprintf('%.3g', u(x==xc,t==tc))};
    
    h = (b-a)/n;
    k = T/m;
    G = max([abs(1 + d*k - 2*c^2*k/h^2 + abs(v*k/h + 1i*2*c^2*k/h^2))
             abs(1 + d*k - 2*c^2*k/h^2 - abs(v*k/h + 1i*2*c^2*k/h^2))]);
    
    table3(i,4) = {sprintf('%.3g', G)};
end
table3 = array2table(table3);
table3.Properties.VariableNames = {'n', 'm', 'u', 'G'};
disp(table3);
end