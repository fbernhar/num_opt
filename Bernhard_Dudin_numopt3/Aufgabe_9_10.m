%% Aufgabe 9 & 10
% 
disp('Aufgabe 9 & 10')


x = [0,0];
lambda = [1,1,1,1];

constraints = @(x)[-1+x(1)+x(2);-1+x(1)-x(2);-1-x(1)+x(2);-1-x(1)-x(2)];           
hesse = @(x)[2,0;0,12*(x(2)-0.75)^2];
delta_f = @(x)[2*(x(1)-1.5);4*(x(2)-0.75)^3];

A = [1,1;1,-1;-1,1;-1,-1];
b = - constraints(x);
f = delta_f(x);
H = hesse(x);



d = quadprog(H,f,A,b);
v = A'\((-H*d)-f);    

options = optimoptions('fmincon','Algorithm','sqp');
fun = @(x)((x(1)-1.5)^2 + (x(2)-0.75)^4);
sol_x = fmincon(fun,x,A,b,[],[],[],[],[],options);

c = delta_f(d)+(v(1)*A(1,:)+v(2)*A(2,:))';