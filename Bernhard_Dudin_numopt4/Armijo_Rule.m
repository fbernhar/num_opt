function [output] = Armijo_Rule(func,delta_func, x, d, alpha)
%% Armijo Linesearch
%
%disp('Armijo_Rule')
%% Implementation of the Algorithm
beta = 0.5;
sigma = 10.0e-4;

if func(x+alpha*d) <= func(x)+sigma*alpha*(delta_func(x)'*d)
    while (func(x+alpha*d)<func(x)+sigma*alpha*(delta_func(x)'*d))
        alpha=alpha/beta;
    end
    alpha = alpha * beta;
else
    while (func(x+alpha*d)>func(x)+alpha*sigma*(delta_func(x)'*d))
        alpha=alpha*beta;
    end
end
output = alpha;
end

