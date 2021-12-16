function [xk, residuum, jacob, info] = GaussNewton(fun,xk,xdata,ydata, gradTol, maxIter, returnJacob)
% GaussNewton computes the Least-Squares pproximation from the point
% xCurrent with the Gauss-Newton algorithm
% fun MUST return a column, xdata MUST be a column
% shapes must fulfill fun(xk, xdata) = ydata
info.num_iter = 0;
h = 1e-4;

if ~exist('gradTol', 'var')
    gradTol = 1e-5;
end

if ~exist('maxIter', 'var')
    maxIter = 100;
end

if ~exist('returnJacob', 'var')
    returnJacob = 0;
    jacob = [];
end

while 1
    resid = ydata - fun(xk, xdata);
    r_line_xk = compute_r_line(xk, xdata, h, fun);
    grad_f = 2*(r_line_xk'*resid);
    if max(abs(grad_f)) < gradTol
       residuum = resid'*resid;
       info.stopping_reason = 'x fulfills gradient tolerance';
       if returnJacob
           jacob = r_line_xk';
       end
       return
    elseif info.num_iter > maxIter
       residuum = resid'*resid;
       info.stopping_reason = 'maximum number of iterations reached';
       if returnJacob
           jacob = r_line_xk';
       end
       return
    end
        
    if size(xdata(:), 1) > size(xk(:))
        dk = linsolve(r_line_xk'*r_line_xk, grad_f/2);
    else
        dk = linsolve(r_line_xk', resid);
    end
    
    xk = xk + dk(:);
    info.num_iter = info.num_iter + 1;
end

    function [r_line] = compute_r_line(x, xdata, h, f)
        r_line = zeros([size(x,1), size(xdata,1)]);
        for i = (1:size(x,1))
            x_diff = x;
            x_diff(i) = x_diff(i) - h;
            J_xi = @(xdata)(f(x, xdata)-f(x_diff, xdata))./h;
            r_line(i,:) = J_xi(xdata);
            
        end
        r_line = r_line';
    end

end

