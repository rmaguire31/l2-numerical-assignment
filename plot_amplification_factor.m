function plot_amplification_factor(coeffs, method)

if nargin < 1
    coeffs = rand(10,3);
end
if nargin < 2
    method = {'F' 'C' 'B'};
end
if ~isrow(method)
    method = method(:);
end
if ischar(method)
    method = {method};
end
th = linspace(0, pi, 100);


for i = 1:size(coeffs,1)
    R = 2*coeffs(i,1)-1;
    C = 2*coeffs(i,2)-1;
    S = coeffs(i,3);
    
    m = cellfun(@get_method, method');
    G = abs(1 + R + m*(cos(th) - 1) + 1i*C*sin(th));
    
    thsp = repmat([0 pi], length(m), 1);
    Gsp(:,1:2) = abs(1 + R + [m m].*(cos(thsp) - 1));
    
    Gsp(:,3) = nan;
    thsp(:,3) = nan;
    c = (m.*(1 + R - m))./(C^2 - m.^2);
    Gsp(abs(c)<=1,3) = abs(1 + R + m(abs(c)<=1).*(c(abs(c)<=1) - 1) ...
                           + 1i*C*sqrt(1-c(abs(c)<=1).^2));
    thsp(abs(c)<=1,3) = acos(c(abs(c)<=1));
    
    plot(th/pi, G);
    hold('on');
    scatter(thsp(:)/pi, Gsp(:));
    hold('off');
    legend(method);
    if i < size(coeffs,1)
        keyboard;
    end
end
    function c = get_method(m)
        switch m
            case 'F'
                c = 2*S + C;
            case 'C'
                c = 2*S;
            case 'B'
                c = 2*S - C;
        end
    end
end

