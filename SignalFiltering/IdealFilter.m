function [ A_pure ] = ideal_filter(w, w_p)
    A_pure = zeros(1,length(w));
    
    for i = 1:length(w)
        w_i = abs(sign(w(i)).*(mod(abs(w(i)) + pi, 2.*pi) - pi));
        
        if (w_i >= w_p(1) && w_i <= w_p(2))
            A_pure(i) = 1;
        else
            A_pure(i) = 10^(-7);
        end
    end
end

