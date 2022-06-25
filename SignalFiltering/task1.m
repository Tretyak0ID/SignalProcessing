% filter parameters
w_p = [0 0.6*pi]; % bandwidth
M = 5; % M
w_j = pi * ((0:M) + 0.5)/(M + 1);
K_d = (w_j >= w_p(1) & w_j < w_p(2));
 
% calculation of filter coefficients by formula (4.8)
h = zeros(1, M + 1);
for k = 0:M
   for j = 0:M-1
       h(M - k + 1) = h(M - k + 1) + K_d(j + 1) * cos(pi * (k) * (j + 0.5) / (M + 1)); % calculating coefficients
   end
   h(M - k + 1) = h(M - k + 1)/(M + 1);
end
 
w = 0:0.001:2*pi;
 
% calculating the frequency response, see (4.1) and (4.7)
A = zeros(1,length(w));
for i = 0:length(w)-1
    A(i + 1) = h(M + 1);
    for k = 1:M
        A(i + 1) = A(i + 1) + 2 * h(M - k + 1) .* cos(w(i + 1) * k); % see (4.7) 
    end
    A(i + 1) = A(i + 1) * exp(-1i * w(i + 1) * M);  % see (4.1) 
end

%Ploting АЧХ filter
hold on; grid on; xlabel('\omega'); ylabel('|K(\omega)|'); title('АЧХ');
plot(w, abs(A))
plot(w, IdealFilter(w, w_p), 'r')
legend('filter', 'ideal filter')

%Ploting АЧХ filter in db
figure; hold on; grid on; xlabel('\omega'); ylabel('|K(\omega)|'); title('AЧХ в децибелах')
plot(w, db(abs(A)))
plot(w, db(abs(IdealFilter(w, w_p))), 'r')
legend('filter', 'ideal filter')

%Ploting ФЧХ filter

    %Calculate ФЧХ
phi = [];
for i=1:length(w)
    if A(i) < 0
        phi(i) = -M*w(i) - pi;
    else
        phi(i) = -M*w(i) + pi;
    end
    
    if phi(i) >= pi
        while 1
            if phi(i) <= pi
                break;
            end

            phi(i) = phi(i) - 2 * pi;
        end
    end


    if phi(i) < -pi
        while 1
            if phi(i) >= -pi
                break;
            end
            phi(i) = phi(i) + 2 * pi;
        end
    end
end
    %ФЧХ plot
figure; hold on; xlabel('\omega'); ylabel('\phi(\omega)'); title('ФЧХ'); grid on
plot(w, phi);
delta_p = max(abs(abs(A)-IdealFilter(w, w_p)).*(w>w_p(1)).*(w<w_p(2)))
delta_s = max(abs(abs(A)-IdealFilter(w, w_p)).*(w>w_s(1)).*(w<w_s(2)))

