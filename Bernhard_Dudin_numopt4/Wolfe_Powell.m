function alpha = Wolfe_Powell(f,grad,x,d,alpha_init)
%% Wolfe-Powell Schrittweitensteuerung.
% f: Funktion.
% grad: Gradient der Funktion f.
% x: Vektor.
% d: Abstiegsrichtung.
% alpha_init: Startwert alpha.
%disp('Wolfe Powell Verfahren')
%% Implementierung.
alpha=alpha_init;sigma=0.01;gamma=1.5;rho=0.3;bool=0; 

phi_null = f(x);
phi_null_grad = grad(x);
phi_alpha = f(x+alpha*d);
phi_alpha_grad = grad(x+alpha*d);

% Expansion
while (phi_alpha<phi_null+sigma*alpha*phi_null_grad'*d) 
    if (phi_null_grad'*d*rho<=phi_alpha_grad*d') 
        bool = 1;
        break; 
    end
    alpha=gamma * alpha; 
    phi_alpha = f(x+alpha*d);
    phi_alpha_grad = grad(x+alpha*d);
end

if ~bool % Kontraktion
    a = 0; b = alpha; alpha = (a + b) / 2;
    while( ~((phi_alpha <= phi_null + sigma * alpha * phi_null_grad'*d) & ...
        (phi_alpha_grad'*d >= phi_null_grad'*d * rho))  & ((b-a)>1e-14) )
        if (phi_alpha <= phi_null + sigma * alpha * phi_null_grad'*d)
            a = alpha; alpha = (a+b) / 2;
        else
            b = alpha; alpha = (a+b) / 2;
        end
        phi_alpha = f(x+alpha*d);
        phi_alpha_grad = grad(x+alpha*d);
    end
end
end
