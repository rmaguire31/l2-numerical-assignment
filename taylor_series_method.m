function taylor_series_method()
% TAYLOR_SERIES_METHOD Calculate the values for Table 1 and Table 2
%
% COPYRIGHT (C) Russell Maguire 2017

n = [10 20 40 80];

y = {{@(t)-cos(t)  @(t,y)sin(t) @(t,y)cos(t) @(t,y)-sin(t) @(t,y)-cos(t)}
     {@(t)exp(5*t) @(t,y)5*y    @(t,y)25*y   @(t,y)125*y   @(t,y)625*y}
     {@(t)2*t      @(t,y)2      @(t,y)0      @(t,y)0       @(t,y)0}}';
t0 = {0 0 0};
tn = {pi 0.1 1};
y0 = {0 1 0};
odes = struct('y', y, 't0', t0, 'tn', tn, 'y0', y0);

f = {@euler @taylor2 @taylor4};
m = {1, 2, 4};
methods = struct('f', f, 'm', m);

table1 = zeros(length(n), length(methods), length(odes));
table2 = table1;
for i = 1:length(n)
    for j = 1:length(methods)
        for k = 1:length(odes)
            f = odes(k).y(1+(1:methods(j).m));
            h = (odes(k).tn - odes(k).t0)/n(i);
            t0 = odes(k).t0;
            y0 = odes(k).y0;
            
            y_est = methods(j).f(f{:}, h, n(i), t0, y0);
            y_actual = odes(k).y{1}(odes(k).tn);
            
            table1(i, j, k) = y_est(end);
            table2(i, j, k) = abs(y_actual - y_est(end));
        end
    end
end
table1 = array2table(reshape(table1, [], size(table1, 2)));
table1.Properties.VariableNames = {'euler', 'taylor2', 'taylor4'};
table1.Properties.RowNames = {'f(x,y)=sin(t),n=10'
                              'f(x,y)=sin(t),n=20'
                              'f(x,y)=sin(t),n=40'
                              'f(x,y)=sin(t),n=80'
                              'f(x,y)=5y(t),n=10'
                              'f(x,y)=5y(t),n=20'
                              'f(x,y)=5y(t),n=40'
                              'f(x,y)=5y(t),n=80'
                              'f(x,y)=2,n=10'
                              'f(x,y)=2,n=20'
                              'f(x,y)=2,n=40'
                              'f(x,y)=2,n=80'};
table2 = array2table(reshape(table2, [], size(table2, 2)));
table2.Properties.VariableNames = {'err_euler', 'err_taylor2', 'err_taylor4'};
table2.Properties.RowNames = {'f(x,y)=sin(t),n=10'
                              'f(x,y)=sin(t),n=20'
                              'f(x,y)=sin(t),n=40'
                              'f(x,y)=sin(t),n=80'
                              'f(x,y)=5y(t),n=10'
                              'f(x,y)=5y(t),n=20'
                              'f(x,y)=5y(t),n=40'
                              'f(x,y)=5y(t),n=80'
                              'f(x,y)=2,n=10'
                              'f(x,y)=2,n=20'
                              'f(x,y)=2,n=40'
                              'f(x,y)=2,n=80'};
disp(table1);
disp(table2);
end