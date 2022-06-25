function signal = GenerateRandomSignal(count_set, x_min, x_max)
%Генерация случайного равномерно распределенного на отрезке [x_min, x_max]
%сигнала, заданного на отрезке count_set + построение графика.

    signal = unifrnd(x_min,x_max,1,length(count_set));
    
    hold on;grid on;
    plot(count_set,signal);
    title('Random signal');
    xlabel('t'); ylabel('x(t)');
end

