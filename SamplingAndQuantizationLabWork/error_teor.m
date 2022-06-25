function [ e ] = error_teor(t, d, sigma, m) 
syms x 
e = 0; 
for j = 1:(length(t)-1) 
    fx = @(x)((x-d(j)).^2).*(1/(sigma.*sqrt(2.*pi))).*exp(-((x-m).^2)./(2.*sigma.^2)); 
    integ = integral(fx, t(j), t(j+1)); 
    e = e + integ; 
end 
end
