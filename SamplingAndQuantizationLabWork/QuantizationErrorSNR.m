function  error_empiric = QuantizationErrorSNR(original_signal, quantorization_signal, x_min, x_max)
%Анализ ошибки квантования: вычисление SNR ошибки
%равномерного квантования, построение графика

    error_empiric   = zeros(1,length(quantorization_signal(:,1)));
    
    for k = 1:1:length(quantorization_signal(:,1))
        A_signal           = sqrt(sum(original_signal.^2))/length(original_signal);
        A_noise            = sqrt(sum((original_signal-quantorization_signal(k,:)).^2))/length(original_signal);
        error_empiric(k)   = 20*log10(A_signal/A_noise);
        error_theoretic(k) = ((x_max-x_min)/2^k)^2/12;
    end
    
    hold on; grid on;
    plot([1:length(error_empiric)], error_empiric);
    xlabel('bits'); ylabel('error');
    title('SNR error');
end

