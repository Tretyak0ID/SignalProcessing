function output_image = BruteForceGouging(input_image, k)
%Грубое уменьшение частоты дискретизации изображения в n раз и вывод
%полученного изображения.

    output_image = zeros(floor(size(input_image,1)/k), floor(size(input_image,2)/k), 3, 'uint8'); 
    for i=1:size(output_image,1)
        for j=1:size(output_image,2)
            output_image(i, j, :) = input_image(floor(i*k), floor(j*k), :);
        end
    end
    imshow(output_image)
end

