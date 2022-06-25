function quant = OptimalLloydMaxQuantizer(original_signal, t, d, count_set)
%Производит оптимальное квантование Ллойда-Макса с заданным числом битов на
%один отсчет, которое определяется в зависимости от передаваемых параметров
%квантователя
    
    len = length(original_signal);
    quant = zeros(1, len);
    for i=1:len
        v = sum(original_signal(i) > t);
        quant(i) = d(v);
    end
    
    hold on; grid on;
    stairs(count_set, quant)
    title('Lloyd-Maks, bits = ', log2(length(d)))
    xlabel('t'); ylabel('quantization x(t)')
end

