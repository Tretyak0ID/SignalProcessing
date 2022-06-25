function [dot, disc_signal] = SignalSamp(sampling_frequency, t_min, t_max, nu)
%Функция, производящая дискретизацию сигнала, при этом дискретизированный
%сигнал должен быть реализован в виде m-функции GenerateSignal

    t_step      = 1/sampling_frequency;   
    dot         = [t_min:t_step:t_max];
    disc_signal = GenerateSignal(t_min, t_max, t_step, nu);
    
    hold on; grid on;
    plot(dot, disc_signal, '.m');
    title('Time domain, sampling frequency ', sampling_frequency);
    legend('original signal', 'sampled signal', 'sampled dot')
end

