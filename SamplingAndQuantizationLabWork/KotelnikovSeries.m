function series = KotelnikovSeries(set_of_samples, t_min, t_max, original_step, sampling_frequency)
%Функция восстанавливает сигнал по набору его отсчетов при помощи теоремы
%Котельникова (предполагает предварительное построение графиков исходного 
%сигнала до дискретизации и набора отсчетов внешним образом)

    syms t

    t_step = 1/sampling_frequency;  
    kt     = [t_min: t_step: t_max];
    KotelnikovSer = matlabFunction(sum(set_of_samples.*(sin(2*pi*sampling_frequency*(t-kt)))./((2*pi*sampling_frequency*(t-kt)))));
    
    t      = [t_min: original_step :t_max];
    series = KotelnikovSer(t);
    
    hold on; grid on;
    plot(t,series);
    title('Time domain')
    xlabel('t'); ylabel('f(t)');
    legend('original signal', 'points of samples', 'alternative signal')
end

