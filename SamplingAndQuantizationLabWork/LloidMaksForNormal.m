function [t_new, d_new] = LloidMaksForNormal(t, d, m, sigma)
%Вычисляет параметры квантователя Ллоида-Макса для нормального
%распределения с произвольными параметрами
    t_new = zeros(1,length(t));
    d_new = zeros(1,length(d));

    for k = 1:length(t)
        t_new(k) = t(k)*sigma + m;
    end
    
    for k = 1:length(d)
        d_new(k) = d(k)*sigma + m;
    end
end

