function frequency_domain = FourierTransformSamp(time_domain, sampling_frequency, F_max)
%Производит преобразование Фурье для дискретного сигнала, заданного набором
%отсчетов во временной области.

    t_step           = 1/sampling_frequency;
    nu               = [-F_max-10: 2*F_max/1000 :F_max+10];
    frequency_domain = zeros(1,length(nu));
    
    for k = 1:length(time_domain)
        frequency_domain = frequency_domain + time_domain(k)*exp(-2*pi.*nu*(k-1)*t_step*j)*t_step;
    end
    
    %АЧХ
    hold on; grid on;
    plot(nu,abs(frequency_domain))
    xlabel('\nu'); ylabel('|S(\nu)|');
    title('Frequency domain')
end

