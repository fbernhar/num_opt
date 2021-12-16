function [f] = regressionFunc2(t)
    syms x1 x2 x3;
    f = x1*exp(-(x2^2 + x3^2)*t).*sinh((x3^2)*t)/(x3^2);
end

