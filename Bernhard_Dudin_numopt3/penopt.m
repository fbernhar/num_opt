function res = penopt(f,h,x0,r,iter, rMax, hTol)
    % f - Funktion zu minimieren
    % h - Gleichheitsnebenbedingungen
    % x0 - startwert
    % r - penalty-parameter
    % iter - maximale Anzahl an Iterationen
    % rMax - obere Grenze für r
    % hTol - erlaubte toleranz zwischen h und 0
    
    %(i) init
    k = 1;
    x(k,:) = x0;
    no_progress_count = 0;
    hxk_old = inf;
    [P,gradient] = getFunAndGrad(f,h);
    while k < iter
        %(ii) Lösung des Minimierungsproblems
        xk = BFGS_Pen(@(x,r) evaluate(P,gradient,x,r),r,x0,100);
        xkk = xk.x(end,:);
        x(k+1,:) = xkk;
        %(iii) falls h(x)=0 stopp
        [~,n] = size(h);
        hi = 0;
        for i = 1:n
            hi = hi + abs(h(i));
        end
        mhx = matlabFunction(hi);
        hxk_check = num2cell(xkk);
        hxk_new = mhx(hxk_check{:});
        if hxk_new < hTol 
            disp("h Toleranz erfüllt")
            res.x = x;
            break;
        end
        if hxk_new < hxk_old
            hxk_old = hxk_new;
            no_progress_count = 0;
        else
            no_progress_count = no_progress_count + 1;
        end

        if no_progress_count >= 5
            disp("Keine Verbesserung nach 5 iterationen")
            break;
        end

        %r erhöhen
        if r < rMax
            r = r * 10;
        else
            disp("maximale Bestrafungsstärke erreicht")
            res.x = x;
            break;
        end
        k = k+1;
    end
    res.x = x;
end

function [mf,mg] = getFunAndGrad(f,h)
    syms r
    p = f;
    for i = 1:size(h)
        p=p+r/2*abs(h(i))^2;
    end
    vars = symvar(p);
    diff_vars = vars(2:end);
    g = jacobian(p,diff_vars);
    mf = matlabFunction(p);
    mg = matlabFunction(g);
end

% f und gradient evaluieren
function [f,g] = evaluate(mf,mg,x,r)
    xn = num2cell(x);
    f = mf(r,xn{:});
    g = mg(r,xn{:});
end