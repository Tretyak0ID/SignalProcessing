function output_signal = UniformQuantization(original_signal, count_set, bits, x_min, x_max)
%Функция равноерного квантования дискретного сигнала и построение графика 
%квантованного сигнала

        N                   = 2^bits;
        q_step              = (x_max-x_min)/N;
        output_signal       = zeros(1,length(count_set));
        
        for i=1:length(original_signal)
            l = floor(original_signal(i)./q_step);
            if ((mod(original_signal(i), q_step) == 0) && (l > 0))
                l = l - 1;
            end
            output_signal(i) = q_step.*(2.*l + 1)./2;
        end
        
        hold on; grid on;
        stairs(count_set, output_signal)
        title('Uniform quantization, bits = ', bits)
        xlabel('t'); ylabel('quantization x(t)')
end

