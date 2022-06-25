w_p = [0 0.6*pi] ; w_s = [0.8*pi pi];  % bandwidth and stopband
sigma_p = 0.0125; sigma_s = 0.015; %delta
delta_p = zeros(1,36); delta_s = zeros(1,36);

for M = 5:40
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

    % calculating the frequency response is similar to the task1
    C = zeros(1,length(w));
    for i = 0:length(w)-1
        C(i + 1) = h(M + 1);
        for k = 1:M
            C(i + 1) = C(i + 1) + 2 * h(M - k + 1) .* cos(w(i + 1) * k); % see (4.7) 
        end
        C(i + 1) = C(i + 1) * exp(-1i * w(i + 1) * M);  % see (4.1) 
    end   

    delta_p(M - 4) = max(abs(abs(C)-IdealFilter(w, w_p)).*(w>w_p(1)).*(w<w_p(2)));
    delta_s(M - 4) = max(abs(abs(C)-IdealFilter(w, w_p)).*(w>w_s(1)).*(w<w_s(2)));
    
    if (delta_p(M-4)<sigma_p)&(delta_s(M-4)<sigma_s)
        M
        break;
    end     
end

sigma_s
delta_s(1:M-3)
sigma_p
delta_p(1:M-3)

%Ploting АЧХ filter
hold on; grid on; xlabel('\omega'); ylabel('|K(\omega)|'); title('АЧХ');
plot(w, abs(C))
plot(w, IdealFilter(w, w_p), 'r')
legend('optimal filter', 'ideal filter')

%Ploting АЧХ filter in db
figure; hold on; grid on; xlabel('\omega'); ylabel('|K(\omega)|'); title('AЧХ в децибелах')
plot(w, db(abs(C)))
plot(w, db(abs(IdealFilter(w, w_p))), 'r')
legend('optimal filter', 'ideal filter')