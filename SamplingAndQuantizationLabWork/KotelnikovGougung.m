function output_image = KotelnikovGougung(input_image,k)

    fd = input_image;
    [M, N] = size(fd);
    Mk = floor(M / k);
    Nk = floor(N / k);
    ff = fd(1:k:(Mk * k), 1:k:(Nk * k));
    ColumnInd = [1:max(Mk, Nk)];
    for j = 1:max(M,N)
        SincArray(j, ColumnInd) = sinc(j / k - ColumnInd);
    end
    F = SincArray(1:M, 1:Mk)*ff*SincArray(1:N, 1:Nk)';
    
    output_image = F;
    imshow(output_image);
end

