%% Project 04
% Florian Bernhard 2825717

%% Aufgabe 01 & 02

x0 = [0, 0];
x1 = [-1.2, 1];
x2 = rand(1,2);

%%Rosenbrock Function.
func_rosenbrock = @(x) [(1-x(1)).^2 + 100*(x(2)-x(1).^2).^2];
delta_rosenbrock = @(x) [2*(1-x(1))+200*(x(2)-x(1).^2)*-2*x(1);200*(x(2)-x(1).^2)];
hesse_rosenbrock = @(x) [-2-400*x(2)+1200*x(1).^2, -400*x(1); -400*x(1), 200];

sol_rosenbrock01 = global_newton_method(func_rosenbrock, delta_rosenbrock, hesse_rosenbrock, x0);
sol_rosenbrock02 = global_newton_method(func_rosenbrock, delta_rosenbrock, hesse_rosenbrock, x1);
sol_rosenbrock03 = global_newton_method(func_rosenbrock, delta_rosenbrock, hesse_rosenbrock, x2);

sol_fminunc01 = fminunc(func_rosenbrock, x0);
sol_fminunc02 = fminunc(func_rosenbrock, x1);

%%Himmelblau Function.
func_himmelblau = @(x) [(x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 -7).^2];
delta_himmelblau = @(x) [2*(x(1).^2 +  x(2) - 11)*2*x(1) + 2*(x(1)+x(2).^2-7);2*(x(1).^2 +  x(2) - 11) + 2*(x(1)+x(2).^2-7)*2*x(2)];
hesse_himmelblau = @(x) [12*x(1)^2-4*x(2)-22+2, -4*x(1)+ 4*x(2); 4*x(1)+4*x(2), 2+ 4*x(1)+12*x(2)^2-28];

sol_himmelblau01 = global_newton_method(func_himmelblau, delta_himmelblau, hesse_himmelblau, x0);
sol_himmelblau02 = global_newton_method(func_himmelblau, delta_himmelblau, hesse_himmelblau, x1);
sol_himmelblau03 = global_newton_method(func_himmelblau, delta_himmelblau, hesse_himmelblau, x2);

sol_fminunc03 = fminunc(func_himmelblau, x0);
sol_fminunc04 = fminunc(func_himmelblau, x1);

%% Aufgabe 06
Ab = [1 1 0 0;0 1 1 0;1 0 0 0; 0 1 0 1];
An = [1 1 0;3 0 0;0 0 1;0 0 0];
Gn_b = (Ab)\An
x = (Ab)\b'
c_b=[-2 -4 0 0];
func_value = c_b*x
cn = [-3 0 0];
Se = c_b*Gn_b-cn
pivot = x./Gn_b
pivot = pivot(:,end)
%% Aufgabe 07

options = optimoptions('linprog','Algorithm','dual-simplex','Display','iter')
A = [1 1 1; 0 3 1;1 0 0;0 0 1];
f = [-2 -3 -4];
b = [4 6 2 3];
linprog(f,A,b,[],[],[],[],options)

Ab = [1 1 0 0;0 1 1 0;1 0 0 1;0 1 0 0];
An = [1 1 0; 3 0 0;0 0 0;0 0 1];
  
Gn_b = (Ab)\An
x = (Ab)\b'
c_b=[-2 -4 0 0]
func_value = c_b*x
cn = [-3 0 0]
Se = c_b*Gn_b-cn
pivot = x./Gn_b
pivot = pivot(:,end)



%% Aufgabe 08

Q = [2 -2; -2 4];
q = [-2 -6];

U = [1/2 1/2; -1 2];
r = [1 2];

sol_active_set = active_set_quadprog(Q,q,[],[],U,r,[])

options = optimoptions(@quadprog,'Display','iter');
sol_quadprog = quadprog(Q,q,U,r,[],[],[],[],[],options)


%% Aufgabe 12
Q = [4 0 0 0;
     0 7/2 0 0;
     0 0 3/2 0
     0 0 0 0];
q = [4500 4000 3500 3000];

U = [1 1 1 0;
     -1 -1 -1 0];
r = [10500 -8500];
G = [1 1 1 1];
b = [10000];

x = [3000 4000 2000 1000];
active_set_quadprog(Q,q,G,b,U,r,x)
options = optimoptions('quadprog','Algorithm','active-set');
quadprog(Q,q,U,r,G,b,[],[],x,options)