function output_signal = GenerateSignal(t_min,t_max,t_step,nu)
%Функция генерации сигнала представляющего собой сумму синусов с набором
%частот nu (массив частот)

    t  = [t_min:t_step:t_max];
    output_signal = sin(2*pi*nu(1).*t) + sin(2*pi*nu(2).*t);
    
    hold on; grid on;
    plot(t, output_signal)
    title('Sum of sines'); xlabel('t'); ylabel('signal');
end

