%Task-2-filter-realized----------------------------------------------------
w_p = [0 0.6*pi] ; w_s = [0.8*pi pi];  % bandwidth and stopband
M = 5; % M
 
w_j = pi * ((0:M) + 0.5)/(M + 1);
flag_p = (w_j > w_p(1)) & (w_j <= w_p(2));
flag_s = (w_j >= w_s(1) & w_j < w_s(2));
N=length(w_j)-sum(flag_p)-sum(flag_s)-1;
x = w_p(2):(w_s(1)-w_p(2))/N:w_s(1); % initial value(s) of the frequency response in the transition band
x = fminsearch(@(x) syntez(x, w_p, w_s, M), x);
[Er, h] = syntez(x, w_p, w_s, M);
 
for k = 1:M
    h(M + k + 1) = h(M - k + 1);
end
%--------------------------------------------------------------------------
T_max = 70;

figure
grid on; hold on
plot(1:T_max, filter(h, 1, sin(0.9*pi*[1:T_max])))
plot(1:T_max, sin(0.9*pi*[1:T_max]))
legend('filtered sin(0.9\pin)', 'unfiltered sin(0.9\pin)')

figure
grid on;hold on
plot(1:T_max, filter(h, 1, sin(0.9*pi*[1:T_max])))
plot(1:T_max, sin(0.9*pi*([1:T_max]-M)))
legend('filtered sin(0.9\pin)', 'unfiltered sin(0.9 \pi(n-M))')

figure
grid on; hold on
plot(1:T_max, filter(h, 1, sin(0.5*pi*[1:T_max])))
plot(1:T_max, sin(0.5*pi*[1:T_max]))
legend('filtered sin(0.5\pin)', 'unfiltered sin(0.5\pin)')

figure
grid on; hold on
plot(1:T_max, filter(h, 1, sin(0.5*pi*[1:T_max])))
plot(1:T_max, sin(0.5*pi*([1:T_max]-M)))
legend('filtered sin(0.5\pin)', 'unfiltered sin(0.5\pi(n-M))')