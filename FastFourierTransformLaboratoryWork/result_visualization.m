%Анализ эффективности БПФ и ДПФ
DFT = [1e-06, 1.3e-05, 4e-05, 0.00015, 0.000872, 0.002446, 0.008645, 0.030223, 0.111084, 0.439266, 1.71226, 6.76674, 27.0718];
FFT = [1e-06, 1.1e-05, 1.8e-05, 3.7e-05, 8.9e-05, 0.000202, 0.000539, 0.001002, 0.001976, 0.004069, 0.008045, 0.017287, 0.035784];
N   = zeros(1,13);

for i=1:13
    N(i) = 2^(i+1);
end

hold on; grid on;
plot(N, DFT, 'm')
plot(N, FFT, 'b')
xlabel('N'); ylabel('t(N)'); legend('time DFT', 'time FFT')

figure
hold on; grid on;
plot(N, log2(DFT), 'm')
plot(N, log2(FFT), 'b')
xlabel('N'); ylabel('log2(t(N))'); legend('time DFT', 'time FFT')

%Анализ эффективности свертки
convol     = [1e-06, 6e-06, 2e-05, 7.3e-05, 0.000238, 0.001114, 0.004074, 0.014295, 0.051315, 0.197949, 0.81008, 3.20245, 12.618];
fft_convol = [1.2e-05, 6.1e-05, 0.000143, 0.000311, 0.000615, 0.001314, 0.002962, 0.00568, 0.011798, 0.026213, 0.057111, 0.116641, 0.247899];

figure
hold on; grid on;
plot(N, convol, 'm')
plot(N, fft_convol, 'b')
xlabel('N'); ylabel('t(N)'); legend('time convolution', 'time fft convolution')

figure
hold on; grid on;
plot(N, log2(convol), 'm')
plot(N, log2(fft_convol), 'b')
xlabel('N'); ylabel('log2(t(N))'); legend('time convolution', 'time fft convolution')

convol_fix     = [1.8e-05, 3.4e-05, 6.6e-05, 0.000133, 0.000283, 0.000648, 0.001627, 0.003891, 0.011765, 0.036265, 0.132484, 0.496657, 1.94634];
fft_convol_fix = [0.000646, 0.000731, 0.000702, 0.000709, 0.000686, 0.001334, 0.002931, 0.006212, 0.012286, 0.025789, 0.054451, 0.122614, 0.269964];

figure
hold on; grid on;
plot(N, convol_fix, 'm')
plot(N, fft_convol_fix, 'b')
xlabel('N'); ylabel('t(N)'); legend('time convolution', 'time fft convolution')

figure
hold on; grid on;
plot(N, log2(convol_fix), 'm')
plot(N, log2(fft_convol_fix), 'b')
xlabel('N'); ylabel('log2(t(N))'); legend('time convolution', 'time fft convolution')





