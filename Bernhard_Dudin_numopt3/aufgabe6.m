f = @(x) x(1);
x0 = [0.5, 0.5];
%optimieren mit fmincon
fmin_res = fmincon(f,x0,[],[],[],[],[],[],@cons1);
disp("Kleiner als - Variante")
disp("x=("+num2str(fmin_res)+")")
[c,~]=cons1(fmin_res);
disp("g1(x)="+num2str(c(1)))
disp("g2(x)="+num2str(c(2)))
disp("f(x)="+f(fmin_res))

x0 = [1, 0];
%optimieren mit fmincon
fmin_res = fmincon(f,x0,[],[],[],[],[],[],@cons2);
disp("Ist gleich - Variante")
disp("x=("+num2str(fmin_res)+")")
[~,ceq]=cons2(fmin_res);
disp("g1(x)="+num2str(ceq(1)))
disp("g2(x)="+num2str(ceq(2)))
disp("f(x)="+f(fmin_res))





function [c,ceq] = cons1(x)
    c(1) = x(1)^2+x(2)^2-1;
    c(2) = x(1)+x(2)-1;
    ceq = [];
end
function [c,ceq] = cons2(x)
    ceq(1) = x(1)^2+x(2)^2-1;
    ceq(2) = x(1)+x(2)-1;
    c = [];
end