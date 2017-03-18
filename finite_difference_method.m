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
f = @(x)4*x + 4*x.^2;

n = [4, 16 32];
m = [10 100 200 400 800];

% Cartesian product of n and m.
[N, M] = meshgrid(n, m);
nm = [N(:) M(:)];

table3 = zeros(size(nm,1), 3);
for i = 1:size(nm)
    n = nm(i,1);
    m = nm(i,2);
    [x, t, u] = conv_diff(a, b, n, T, m, c, v, d, f);
    
    surf(t, x, u);
    
    table3(i,1) = n;
    table3(i,2) = m;
    table3(i,3) = u(x==xc,t==tc);
end
table3 = array2table(table3);
table3.Properties.VariableNames = {'n', 'm', 'u'};
disp(table3);
end