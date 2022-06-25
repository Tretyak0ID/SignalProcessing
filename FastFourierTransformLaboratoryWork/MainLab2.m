n = 4; %Размерность векторов для преобразований Фурье и сверток

%ГЕНЕРАЦИЯ-СИГНАЛОВ-ДЛЯ-ФУРЬЕ-ПРКОБРАЗОВАНИЯ------------
t = linspace(-5,5,2^n);
st = -1000; en = 1000;
x = randi([st en], 1, length(t))/100;
y = randi([st en], 1, length(t))/100;
%сохранение в файл
fid = fopen('data/real.txt', 'wt');
fprintf(fid, '%f ', x);
fclose(fid);
type('data/real.txt');

fid = fopen('data/imag.txt', 'wt');
fprintf(fid, '%f ', y);
fclose(fid);
type('data/imag.txt');

%МАТЛАБОВСКОЕ-ФУРЬЕ-ПРЕОБРАЗОВАНИЕ---------------------
vect1     = x + j.*y;
fft_vect1 = fft(vect1)./sqrt(length(vect1));

fid = fopen('data/matlab_fft_real.txt', 'wt');
fprintf(fid, '%f ', real(fft_vect1));
fclose(fid);
type('data/matlab_fft_real.txt');
fid = fopen('data/matlab_fft_imag.txt', 'wt');
fprintf(fid, '%f ', imag(fft_vect1));
fclose(fid);
type('data/matlab_fft_real.txt');


%ГЕНЕРАЦИЯ-СИГНАЛОВ-ДЛЯ-СВЕРТОК-------------------------
t = linspace(-5,5,2^n);
t2 = linspace(-5,5,2^6);
st = -1000; en = 1000;
x1 = randi([st en], 1, length(t))/100;
y1 = randi([st en], 1, length(t))/100;
x2 = randi([st en], 1, length(t2))/100;
y2 = randi([st en], 1, length(t2))/100;
%сохранение в файл
fid = fopen('data/real_convolution_x.txt', 'wt');
fprintf(fid, '%f ', x1);
fclose(fid);
type('data/real_convolution_x.txt');

fid = fopen('data/imag_convolution_x.txt', 'wt');
fprintf(fid, '%f ', y1);
fclose(fid);
type('data/imag_convolution_x.txt');

fid = fopen('data/real_convolution_y.txt', 'wt');
fprintf(fid, '%f ', x2);
fclose(fid);
type('data/real_convolution_y.txt');

fid = fopen('data/imag_convolution_y.txt', 'wt');
fprintf(fid, '%f ', y2);
fclose(fid);
type('data/imag_convolution_y.txt');

%МАТЛАБОВСКИЕ-СВЕРТКИ---------------------------------
vect_conv_1 = x1 + j.*y1;
vect_conv_2 = x2 + j.*y2;
vect_conv   = conv(vect_conv_1, vect_conv_2);

fid = fopen('data/matlab_conv_real.txt', 'wt');
fprintf(fid, '%f ', real(vect_conv));
fclose(fid);
type('data/matlab_conv_real.txt');
fid = fopen('data/matlab_conv_imag.txt', 'wt');
fprintf(fid, '%f ', imag(vect_conv));
fclose(fid);
type('data/matlab_conv_real.txt');
