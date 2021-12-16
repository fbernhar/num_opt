function [mu] = Bisektion(func,lower_bound,upper_bound)
disp('Bisektions-Verfahren')

% Auswertung der Funktion an den Grenzen des Intervalls.
aVal = feval(func,lower_bound);
bVal = feval(func,upper_bound);

% Maximale Iterationsschritte festlegen.
i = 0;
it_max = 50;


disp('iter     lower_bound        x0           upper_bound       f(x0)');
while(i < it_max)
    i = i + 1;
    
    %Mittelpunkt des Intervalls bestimmen.
    mu = (lower_bound + upper_bound)/2;
    lambda = feval(func,mu);
    
    % Konsolenausgabe.
    if(length(mu) == 1)
        fprintf('%2i \t %f \t %f \t %f \t %f \n', i-1, lower_bound, mu, upper_bound, lambda);
    elseif(length(mu) == 2)
        fprintf('%2i \t %f %f \t %f %f \t %f %f \t %f \n', i-1, lower_bound, mu, upper_bound, lambda);
    end
    
    %Intervallgrenzen änderungen.
    if aVal > lambda %&& aVal > bVal
        lower_bound = mu;
        aVal = lambda;
    else
        upper_bound = mu;
        bVal = lambda;
    end;
    
end
%Letze Intervalländerung mit anschliender Konsolenausgabe.
mu = (lower_bound + upper_bound)/2;
lambda = feval(func,mu);
if(length(mu) == 1)
    fprintf('\n x = %f produces f(x) = %f \n with %i iterations\n', mu, lambda, i);
elseif(length(mu) == 2)
    fprintf('\n x = %f %f produces f(x) = %f \n with %i iterations\n', mu, lambda, i);
end;
end


