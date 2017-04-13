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
    
    R = d*k;
    C = v*k/h;
    S = c^2*k/h^2;
    
    Gsp = abs(1 + R + (2*S+C)*([-1 1] - 1));
    c3 = ((2*S+C)*(1 + R - (2*S+C)))/(C^2 - (2*S+C)^2);
    if (abs(c3) <= 1)
        Gsp(end+1) = abs(1 + R + (2*S+C)*(c3 - 1) + 1i*C*sqrt(1-c3^2));
    end
    G = max(Gsp);
    
    table3(i,4) = {sprintf('%.3g', G)};
end
table3 = array2table(table3);
table3.Properties.VariableNames = {'n', 'm', 'u', 'G'};
disp(table3);
end