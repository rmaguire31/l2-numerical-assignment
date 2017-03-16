function [x,t,u] = conv_diff(a,b,n,T,m,c,v,d,f)
%%

h = (a-b)/n;
k = T/m;

x = a:h:b;
t = 0:k:T;

diags = repmat(k/h * [c^2+v, -2*c^2-v+d*h, c^2], [min(m, n) 1]);
Dx = spdiags(diags, -1:1, m, n);
full(Dx)

end