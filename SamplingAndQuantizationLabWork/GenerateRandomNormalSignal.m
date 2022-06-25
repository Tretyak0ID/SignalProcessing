function output_signal = GenerateRandomNormalSignal(m, sigma, count_set)
%Генерация случайного дискретного равномерно распределенного случайного
%сигнала.

    output_signal = m + sqrt(sigma).*randn(1, length(count_set));
    
    hold on; grid on;
    plot(count_set, output_signal);
    xlabel('t'); ylabel('x(t)');
    title('normal signal');
end

