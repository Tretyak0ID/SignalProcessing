%Task-2-filter-realized----------------------------------------------------
w_p = [0 0.6*pi] ; w_s = [0.8*pi pi];  % bandwidth and stopband
M = 12; % M
 
w_j = pi * ((0:M) + 0.5)/(M + 1);
flag_p = (w_j > w_p(1)) & (w_j <= w_p(2));
flag_s = (w_j >= w_s(1) & w_j < w_s(2));
N=length(w_j)-sum(flag_p)-sum(flag_s)-1;
x = w_p(2):(w_s(1)-w_p(2))/N:w_s(1); % initial value(s) of the frequency response in the transition band
x = fminsearch(@(x) syntez(x, w_p, w_s, M), x);
[Er, h] = syntez(x, w_p, w_s, M);
 
for k = 1:M
    h(M + k + 1) = h(M - k + 1);
end
%--------------------------------------------------------------------------
im = im2double(imread('var1.png'));

image   = zeros(450+M,297+M);
image(1:450,1:297) = im;

%Filtered image mass
filt_rows = zeros(size(image));
filt_columns = zeros(size(filt_rows));
diff = zeros(size(filt_columns));

%Filtering
for i = 1:size(image,1)
    filt_rows(i,:) = filter(h, 1, image(i,:));
end

for i = 1:size(image,2)
    filt_columns(:, i) = filter(h, 1, filt_rows(:, i));
end

%Show image
figure;
imshow(image);

filt_columns(1:450,1:297) = filt_columns(M+1:450+M,M+1:297+M)
filt_columns(451:end,298:end) = 0;
diff = filt_columns - image;

figure;
imshow(filt_columns);

figure;
imshow(diff); 