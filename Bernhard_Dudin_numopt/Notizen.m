% 1.Funktion: Himmelblau Funktion
f = @(x)( (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2 );
a = [2, 1];
b = [4, 3];

tic
x_bisektion = Bisektion(f,a,b);
t_bisektion = toc;


tic
x0 = [2,4];
x_mutation = Mutation(f,x0);
t_mutation = toc;


tic
x_minsearch = fminsearch(f,x0);
t_minsearch = toc;


% 2.Funktion: Bazaraa Shetty Funktion
g = @(x)(100*(x(1)-2)^4 + (x(1) - 2*x(2))^2);
a = [2, 1];
b = [4, 3];

tic
x_bisektion = Bisektion(g,a,b);
t_bisektion = toc;

tic
x0 = [4,2];
x_mutation = Mutation(g,x0);
t_mutation = toc;

tic
x_minsearch = fminsearch(g,x0);
t_minsearch = toc;

% 3.Funktion
h = @(x)(exp(-x)+0.5*x^2);
a = 0;
b = 1;

tic
x_bisektion = Bisektion(h,a,b);
t_bisektion = toc;

tic
x0 = 1;
x_mutation = Mutation(h,x0);
t_mutation = toc;

tic
x_minsearch = fminsearch(h,x0);
t_minsearch = toc;

% Aufgabe 5
tic
% in fminsearch2 TolF entspricht der maximal erlaubten Standartabweichung
x_minsearch = fminsearch2(h,x0, optimset('Display','iter','TolF',1e-7));
t_minsearch = toc;

% Sehr flache Parabel als ein Beispiel einer Zielfunktion, fuer die das
% Nelder-Mead Kriterium keine gute ergibnisse liefert
flat_parabola = @(x) (1e-10)*x.*x;
x0 = 3;
disp(newline)
disp('True solution is 0')
res1 = fminsearch2(flat_parabola,x0, optimset('TolF',1e-6));
res2 = fminsearch(flat_parabola,x0, optimset('TolF',1e-6));
disp(['result with fminsearch2, tolerance 1e-6 :', num2str(res1)]);
disp(['result with fminsearch, tolerance 1e-6 :', num2str(res2)]);
res1 = fminsearch2(flat_parabola,x0, optimset('TolF',1e-15));
res2 = fminsearch(flat_parabola,x0, optimset('TolF',1e-15));
disp(['result with fminsearch2, tolerance 1e-15 :', num2str(res1)]);
disp(['result with fminsearch, tolerance 1e-15 :', num2str(res2)]);

