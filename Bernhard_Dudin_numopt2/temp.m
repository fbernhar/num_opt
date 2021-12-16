J2 = @(x,h,F)@(t)(F(repmat(x(:),size(x(:)'))+diag(h), t)-F(repmat(x(:),size(x(:)')),t))./h';
f = @(x,t) [t.*x(1)^2; t.*x(2)^3];
t = (1:10)
h = 1e-3
x = [1;1]
a = J2(x,h,f)
ret = f(x,t)
a(t)