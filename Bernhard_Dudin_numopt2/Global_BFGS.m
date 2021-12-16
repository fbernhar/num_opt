function [x, struct_inf] = Global_BFGS(f, grad, x, tol)
%% Implementierung  des Globalen Quasi-Newton-Verfahrens.
% f:     Ist die übergebene Funktion in Form eines Function Handles.
% grad:  Ist der Gradient der Funktion in Form eines Function Handles.
% x:     Ist der Startvektor.
% tol:   Ist die Toleranz für das Abbruchkriterium der Implementierung.
disp('Globales Quasi-Newton-Verfahren')
%%
k = 0;          %Iterationscounter.
it_inf = [];    %Iteration - Container.
alpha = 1;      %Alpha Initialiserung für die Wolfe Powell Schrittweite.
struct_inf = struct(); % Datenkapselung für die Ausgabe.

f_inf = [];            %Funktionswerte - Container.
grad_inf = [];         %Gradientenwerte - Container.
A = eye(length(x));    %Approximation Hesse Matrix.
A_inf = [];            %Hesse Matrix Approximation - Container.
err = 2 * tol;         %Fehler initialsierung.      
err_inf = [];          %Fehler - Container.

while err > tol
    d = - A * grad(x);
    
    if ((k > 0)&&(y'*s < 0 ))             %(y'*s>0)Notwendige Bedingung an die Approximation von A.
        d = -grad(x);
        %A = (y'*s)/(y'*y)*eye(length(x)); %1. Möglichkeit A neu zu berechnen.
        A = eye(length(x));                %2. Möglichkeit A zu berechnen.
                                           %3. Möglichkeit, man behält das A bei.
    end
    alpha = Wolfe_Powell(f, grad, x, d, alpha);
    x_prev = x;
    x = x_prev + alpha * d;
    
    s = x - x_prev;
    y = grad(x) - grad(x_prev);
    update_term = (((s-A*y)*(s')+s*(s-A*y)')/((y')*s)) ...
        -((((s-A*y)')*y*s*(s'))/((y')*s)^2);
    A = A + update_term;
    err = norm(grad(x));
    k  =  k  + 1; 
    
    % Daten aus den Iterationen in den Containern sammeln.
    it_inf = [it_inf,k];  %#ok<*AGROW>
    f_inf = [f_inf,f(x)];
    grad_inf = [grad_inf,grad(x)];
    A_inf = [A_inf,A];
    err_inf = [err_inf,err];
    
end
% Structfields erzeugen und Informationen der einzelnen Iterationen
% abspeichern.
struct_inf.it_inf = it_inf;
struct_inf.f_inf = f_inf;
struct_inf.grad_inf = grad_inf;
struct_inf.A_inf = A_inf;
struct_inf.err_inf = err_inf;
end

