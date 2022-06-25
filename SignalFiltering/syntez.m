% objective functions
function [Er, h] = syntez(x_j, w_p, w_s, M)
 
w_j = pi * ((0:M) + 0.5)/(M + 1);
flag_p = (w_j > w_p(1)) & (w_j <= w_p(2)); flag_s = (w_j >= w_s(1)) & (w_j < w_s(2));
flag_K_d(flag_s) = 0; flag_K_d(flag_p) = 1;
flag_K_d(~(flag_s|flag_p)) = x_j;
h = zeros(1, M + 1); S = zeros(1, M + 1);
 
for k = 0:M
   for j = 0:M
       S(k + 1) = S(k + 1) + flag_K_d(j+1)*cos(pi * (k) * (j + 0.5) / (M + 1));% it is necessary to add
   end
   h(M - k + 1) = S(k + 1) / (M + 1);
end
 
syms w;
A_w = h(M + 1);
for k = 1:M
    A_w = A_w + 2 * h(M - k + 1) * cos(w * k); 
end
A_w = matlabFunction(A_w);
i = 0;
E = zeros(1,length(w_s(1):0.001:w_s(2)));
for t = w_s(1):0.001:w_s(2)
    i = i + 1;
    E(i) =  - A_w(t);
end
for t = w_p(1):0.001:w_p(2)
    i = i + 1;
    E(i) = 1 - A_w(t);
end
Er = max(abs(E));
end
