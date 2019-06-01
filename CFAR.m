function [I, th] = CFAR(x, R, pfa, Ng)

Nsamp = size(x, 1);
th = zeros(Nsamp, 1);
alpha = R .* (pfa.^(-1/R) - 1);

D = [1:(R-1)/2 - (Ng)/2, (R-1)/2 + 1, (R-1)/2 + Ng/2 + 2 : R];


    for i = 1:Nsamp - R

        %temp = abs(x(i+(R - 1)./2));
        temp = 0;

        for j = D
            temp = temp + abs(x(i + j));
        end
        temp_avg = temp./(R-Ng);

        if i == 1 
            th(1:i + (R - 1)/2) = temp_avg .* alpha;
        elseif i == Nsamp - R
            th(i + (R - 1)/2:Nsamp) = temp_avg .* alpha;
        else
            th(i + (R - 1)/2) = temp_avg .* alpha;
        end
    end
    
    I_i = x(x > th);
    
    I = find(x == max(I_i));
    
end
