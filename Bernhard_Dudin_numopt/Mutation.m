function [x] = Mutation(func,x)
disp('Mutation-Sektion-Verfahren')
% Festlegung einer maximalen Schrittweite über k
k = 1;
k_max = 50;
r = -0.5+(0.5+0.5).*rand(length(x), k_max); % Definieren eines Zufallsvektors

%alpha = 0.7;                    %Konstantes Alpha
alpha = rand(1,k_max);          %Zufällige Alpha
%alpha = linspace(1,0,k_max);    %Kleiner werdenendes Alpha

%count = 0;        Abbruchkriterium Anzahl an Schritten ohne Änderung
%err = 2 * TOL;    Abbruchkriterium Funktionswerte, x Werte
%TOL = 1.0e-10;    
disp('iter           x                 f(x)');
while(k <  k_max ) %&& count < 10 && err > Tol)
    x_k = x  + alpha(k) * r(:,k)';
    
    if func(x_k)<func(x)
        % err = abs(func(x_k)-func(x))  Abbruchkriterium Funktionswerte
        % err = abs(norm(x_k)-norm(x))  Abbruchkriterium x Werte
        x = x_k;
        %count = 0;                     
        %alpha = (1 - 0.1) *  alpha     %Änderung bei Schrittänderung.
    %else
        %count = count + 1;
    end
    if length(x) == 2
        fprintf('%2i \t %f %f \t %f \n', k, x , func(x));
    elseif length(x) == 1
        fprintf('%2i \t %f \t %f \n', k, x , func(x));
    end  
    k = k + 1;
    
end
if(length(x) == 1)
    fprintf('\n x = %f produces f(x) = %f \n with %i iterations\n', x, func(x), k);
elseif(length(x) == 2)
    fprintf('\n x = %f %f produces f(x) = %f \n with %i iterations\n', x, func(x), k);
end;

