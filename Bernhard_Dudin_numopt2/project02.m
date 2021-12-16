%% 2. Kleines Projekt  zur Numerischen Optimierung
% Quasi-Newton-Verfahren & Gaus-Newton-Verfahren
% Florian Bernhard & Maxim Dudin
%% Augabe 4.
% Teil 1.
f = @(x)((x(1).^2 + x(2) -11).^2 + (x(1) + x(2).^2 -7).^2);
grad_f = @(x)([ 2 * (x(1).^2 + x(2) - 11) * 2 * x(1) + 2 * (x(1) + x(2).^2 -7) * 1  ;...
    2 * (x(1).^2 +x (2) - 11)* 1 + 2 * (x(1) + x(2).^2 -7 ) * 2 * x(2) ]);

x0 = rand(2,1);
options = optimset('Display','iter','Algorithm','quasi-newton','HessUpdate','bfgs');
[x_fminunc,fval,exitflag,struct_inf_fminunc] = fminunc(f,x0,options);
[x_global_bfgs, struct_inf_global_bfgs] = Global_BFGS(f,grad_f,x0, 1.0e-10);


%% Aufgabe 4. 
% Teil 2.
N = 10;
syms x [1 N];

x0 = rand(N,1)*10;
g = sum(((1-x(1:N-1)).^2+100*(x(2:N)-x(1:N-1).^2).^2));
grad_g = jacobian(g);

g = matlabFunction(g);
grad_g = matlabFunction(grad_g);

options = optimset('Display','iter','Algorithm','quasi-newton','HessUpdate','bfgs');
[x_fminunc,fval,exitflag,struct_inf_fminunc] = fminunc(@(x) functrans(g,x),x0,options);
[x_global_bfgs, struct_inf_global_bfgs] = Global_BFGS(@(x)functrans(g,x),@(x) functrans(grad_g,x),x0, 1.0e-5);

%% Aufgabe 6
fun = @(x, xdata) x(1).*exp(x(2).*xdata);
x = [1,1]';
xdata = [0,1,2,3]';
ydata = [2.0, 0.7, 0.3, 0.1]';
[xk, residuum, jacob, info] = GaussNewton(fun,x,xdata,ydata, 1e-6, 100, 1);
disp(['stopped after: ', num2str(info.num_iter), ' iterations'])
disp(['reason: ', info.stopping_reason])
disp(['xk: ', mat2str(xk)])
disp(['residuum: ', mat2str(residuum)])
disp(['jacobian: ', mat2str(jacob)])

disp(['solution with  lsqcurvefit:', mat2str(lsqcurvefit(fun,x,xdata,ydata))])

% second test
g = @(x,xdata) x(1).*exp(-(x(2)^2 + x(3)^2).*xdata).*(sinh(x(3)^2.*xdata))./x(3)^2;
xdata = (6:6:180)';
%x = [0.1,0.01,0.01]';
x = [10,0.05, 0.1]';
%x = [50,0.05, 0.1]';

ydata = [24.19, 35.34, 43.43, 42.63, 49.92, 51.53,...
    57.39, 59.56, 55.60, 51.91, 58.27, 62.99, ...
    52.99, 53.83, 59.37, 62.35, 61.84, 61.62, ...
    49.64, 57.81, 54.79, 50.38, 43.85, 45.16, ...
    46.72, 40.68, 35.14, 45.47, 42.40, 55.21]';
[xk, residuum, jacob, info] = GaussNewton(g,x,xdata,ydata, 1e-10, 100, 1);
disp(['stopped after: ', num2str(info.num_iter), ' iterations'])
disp(['reason: ', info.stopping_reason])
disp(['xk: ', mat2str(xk)])
disp(['residuum: ', mat2str(residuum)])
disp(['jacobian: ', mat2str(jacob)])

disp(['solution with  lsqcurvefit:', mat2str(lsqcurvefit(g,x,xdata,ydata))])
y_approx = g(xk,xdata);
figure; hold on
plot(xdata, ydata);
plot(xdata, y_approx);
title('Datensatz und Least Squares mit Gauss-Newton')
hold off
legend('data', 'least squares approximation')


%% Aufgabe 9

t1 = [0; 1; 2; 3];
y1 = [2.0; 0.7; 0.3; 0.1];
t2 = [6; 12; 18; 24; 30; 36; 42; 48; 54; 60; 66; 72; 78; 84; 90; 96; 102; 108; 114; 120; 126; 132; 138; 144; 150; 156; 162; 168; 174; 180];
y2 = [24.19; 35.34; 43.43; 42.63; 49.92; 51.53; 57.39; 59.56; 55.60; 51.91; 58.27; 62.99; 52.99; 53.83; 59.37; 62.35; 61.84; 61.62; 49.64; 57.81; 54.79; 50.38; 43.85; 45.16; 46.72; 40.68; 35.14; 45.47; 42.40; 55.21];


regression_func_1 = regressionFunc1(t1)-y1;
regression_func_1 = regression_func_1'*regression_func_1;
regression_func_2 = regressionFunc2(t2)-y2;
regression_func_2 = regression_func_2'*regression_func_2;

f = matlabFunction(regression_func_1);
f_grad = jacobian(regression_func_1);
f_grad = matlabFunction(f_grad);

g = matlabFunction(regression_func_2);
g_grad = jacobian(regression_func_2);
g_grad = matlabFunction(g_grad);

x0 = [2;-3];
[x_global_bfgs_1, struct_inf_global_bfgs_1] = Global_BFGS(@(x)functrans(f,x),@(x) functrans(f_grad,x),x0, 1.0e-5);

x0 = [2; 0.1; 0.2];
[x_global_bfgs_2, struct_inf_global_bfgs_2] = Global_BFGS(@(x)functrans(g,x),@(x) functrans(g_grad,x),x0, 1.0e-5);

