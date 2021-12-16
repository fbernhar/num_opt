function [output] = global_newton_method(func, delta_func, hesse_func, x0)
%% Globalisiertes Newton-Verfahren
%   func: Die hier betrachtete Funktion 
%   x:    Der Startvektor für die Funktion *func*
%disp('global newton method')
%% Implementierung des Algorithmus
x = x0';
rho = 10.0e-8;
p = 2.1;
beta = 0.5;
sigma = 1.0e-4;
epsilon = 10.0e-12;

iter  = 0;
max_it = 500;
err = 2 * epsilon;
data = struct;
while err > epsilon && iter < max_it
    
    d = hesse_func(x)\-delta_func(x);
    %d = cgs(hesse_func(x),-delta_func(x)); %CG-Method
    
    %Überprüfen ob das System Lösbar ist.
    if ~(rank(hesse_func(x),d) == rank(hesse_func(x)))
        d = -delta_func(x);
    end
    
    %Überprüfung der Bedingung
    if delta_func(x)' * d > - rho*norm(d)^p
        d = -delta_func(x);
    end
    
    %alpha = Armijo_Rule(func,delta_func, x, d, 1);
    alpha = Wolfe_Powell(func,delta_func, x, d, 1);
    x = x + alpha * d;
    
    err = norm(delta_func(x));
    iter = iter +1;
    
end
data.x = x;
data.y = func(x);
data.it = iter;
output = data;
end

