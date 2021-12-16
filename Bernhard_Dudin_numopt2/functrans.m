function [y] = functrans(f,x)
%% Inputargument Transformation Funktion.
%   Wertet den Vektor für die N sym Variablen aus.
%   disp('functrans')
%% Implementation.
x = num2cell(x);
y = f(x{:})';
end

