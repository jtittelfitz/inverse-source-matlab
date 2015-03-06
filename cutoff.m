function [ C ] = cutoff( c,x1,x2,cutoff_value )

    cut = ones(size(c));

    index = abs(x1) > 2 | abs(x2) > 2;
    cut(index) = 0;

    index = x1 > 2 & x1 < 2.5 & abs(x2) <= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) - 2).^2));
    index = x2 > 2 & x2 < 2.5 & abs(x1) <= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x2(index) - 2).^2));
    index = x1 > -2.5 & x1 < -2 & abs(x2) <= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) + 2).^2));
    index = x2 > -2.5 & x2 < -2 & abs(x1) <= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x2(index) + 2).^2));


    %% round corners
    index = x1 > 2 & x1 < 2.5 & x2 > 2 & x2 < 2.5 & (x1 - 2).^2 + (x2 - 2).^2 < 0.5^2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) - 2).^2 - (x2(index) - 2).^2));
    index = x1 > 2 & x1 < 2.5 & x2 > -2.5 & x2 < -2 & (x1 - 2).^2 + (x2 + 2).^2 < 0.5^2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) - 2).^2 - (x2(index) + 2).^2));
    index = x1 > -2.5 & x1 < -2 & x2 > 2 & x2 < 2.5 & (x1 + 2).^2 + (x2 - 2).^2 < 0.5^2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) + 2).^2 - (x2(index) - 2).^2));
    index = x1 > -2.5 & x1 < -2 & x2 > -2.5 & x2 < -2 & (x1 + 2).^2 + (x2 + 2).^2 < 0.5^2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (x1(index) + 2).^2 - (x2(index) + 2).^2));

    %% Square corners:
    %{
    index = x1 > 2 & x1 < 2.5 & x2 >= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (max(x1(index),x2(index)) - 2).^2));
    index = x1 > 2 & x1 < 2.5 & x2 <= -2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (max(x1(index),-x2(index)) - 2).^2));
    index = x1 > -2.5 & x1 < -2 & x2 >= 2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (max(-x1(index),x2(index)) - 2).^2));
    index = x1 > -2.5 & x1 < -2 & x2 <= -2;
    cut(index) = exp(4)*exp(-1./(0.5^2 - (max(-x1(index),-x2(index)) - 2).^2));
    %}

    index = abs(x1) >= 2.5 | abs(x2) >= 2.5;
    cut(index) = 0;
    
    %imagesc(cut)

    C = cut.*c + cutoff_value*(1 - cut);  

    %imagesc(c_33)

end

