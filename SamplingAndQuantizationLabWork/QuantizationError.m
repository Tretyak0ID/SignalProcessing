function  [error_empiric, error_theoretic] = QuantizationError(original_signal, quantorization_signal, x_min, x_max)
%Анализ ошибки квантования: вычисление эмпирической и теоретической ошибок
%равномерного квантования, построение графиков

    error_empiric   = zeros(1,length(quantorization_signal(:,1)));
    error_theoretic = zeros(1,length(quantorization_signal(:,1)));
    
    for k = 1:1:length(quantorization_signal(:,1))
        error_empiric(k)   = sum((original_signal-quantorization_signal(k,:)).^2)./length(original_signal);
        error_theoretic(k) = ((x_max-x_min)/2^k)^2/12;
    end
    
    hold on; grid on;
    plot([1:length(error_empiric)], error_empiric);
    plot([1:length(error_theoretic)], error_theoretic);
    xlabel('bits'); ylabel('error');
    title('Theoretic and empiric error');
    legend('empiric error', 'theoretic error');
end

