%Vorbereitung
%für penopt
syms x1 x2 x3 r
f = 1000-x1^2-2*x2^2-x3^2-x1*x2-x1*x3;
h1 = x1^2+x2^2+x3^2-25;
h2 = 8*x1+14*x2+7*x3-56;
%für fmincon
pf = @(x) 1000-x(1)^2-2*x(2)^2-x(3)^2-x(1)*x(1)-x(1)*x(3);

%startwert
x0 = [2,4,1];

%optimieren mit fmincon
fmin_res = fmincon(pf,x0,[],[],[],[],[],[],@cons);
disp("fmincon")
disp("x=("+num2str(fmin_res)+")")
[~,ceq]=cons(fmin_res);
disp("h1(x)="+num2str(ceq(1)))
disp("h2(x)="+num2str(ceq(2)))
disp("f(x)="+pf(fmin_res))

%optimieren mit penopt
%Nebenbedingungen
neben_b = [h1,h2];

result = penopt(f,neben_b,x0,100,50, 1e15, 1e-6);
res = result.x;
disp("penopt")
disp("x=("+num2str(res(end,:))+")")
disp("h1(x)="+double(subs(h1,{x1,x2,x3},res(end,:))))
disp("h2(x)="+double(subs(h2,{x1,x2,x3},res(end,:))))
disp("f(x)="+double(subs(f,{x1,x2,x3},res(end,:))))
return;

function [c,ceq] = cons(x)
    ceq(1) = x(1)^2+x(2)^2+x(3)^2-25;
    ceq(2) = 8*x(1)+14*x(2)+7*x(3)-56;
    c = [];
end