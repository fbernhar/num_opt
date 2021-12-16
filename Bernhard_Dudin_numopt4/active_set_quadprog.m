function [output] = active_set_quadprog(Q,q,G,b,U,r,x)
%% Active-Set  Methode für quadratische Optimierungsprobleme.
%   Q:  Quad. und gemischte Terme der Funktion
%   q:  alleinstehende Terme der Funktion
%   G: linke Seite, Gleichheitsbedingung
%   b: rechte Seite, Gleichheitsbedingung
%   U: linke Seite, Ungleichheitsbedinung
%   r: rechte Seite, Ungleichheitsbedingung
%   x: Startvektor, falls vorhanden
%disp('active_set_method')
%% Implementation of the Algorithm
equality_constraints = true;
inequality_constraints = true;

W = [];
it = 0;
max_it = 100;
tol = 1.0e-15;
output = struct;

if isempty(G)
    equality_constraints = false;
end
if isempty(U)
    inequality_constraints = false;
end


%Startwertberechnung Phase 1 Simplex-Verfahrens wenn x nicht vorhanden,
if isempty(x)
    UG = [U;G];
    slack = eye(size(UG,1));
    UG_prev = UG;
    UG = [UG, slack];
    for i = (size(U,1)+size(U,2))+1:size(UG,2)
        UG(:,i) = zeros(1,size(UG_prev,1));
    end
    x_var = UG\[r b]';
    x = x_var(1:size(UG_prev,2),:);
    
    for i = size(U,2)+1:size(x_var,1)-size(G,1)
        if x_var(i) == 0
            W(size(W,1)+1) = i-size(UG_prev,2);
        end
    end
    
else % ansonsten alle ieq aktiv.
    x = x';
    for i = 1:size(U,1)
        W(i) = i;
    end
end
% Lösung des einfachen Falls mit Gleichungsnebenbedingungen und ohne
% Ungleichheitsbedingungen, ansonsten
if equality_constraints && ~inequality_constraints
    zero = zeros(size(G,1));
    LinEq = [Q, G';...
             G, zero];
    sol_LinEq = LinEq\[-q b]';
    x = sol_LinEq(1:size(Q,1));
    output.x = [x];
    output.it = [it];
    return;
else % wird hier weiter berechnet.
    while it < max_it
        it = it +1;
        
        delta_func = Q*x+q';
        N = U(W,:);
        N = [G;N];
        zero = zeros(size(N,1));
        
        LinEq = [Q, N';...
                 N, zero];
        if rank(LinEq) < min(size(LinEq)) % Falls die Matrix singulär ist
            % least Squares für eine Näherung verwenden.
            sol = lsqr(LinEq,[-delta_func' zeros(1,size(zero,1))]');
        else
            sol = LinEq\[-delta_func' zeros(1,size(zero,1))]';
        end
        d = sol(1:size(Q,1),:);
        lambda = sol( size(Q,1) + size(G,1) +1 :end ,:);
        [min_lambda,min_index_lambda] = min(lambda);
        
        
        
        %Case 1
        if norm(d) < tol && min_lambda > 0
            output.x = [x];
            output.it = [it];
            return;
        end
        %Case 2
        if norm(d) < tol && min_lambda < 0
            W(min_index_lambda) = [];
            continue;
        end
        %Case 3
        %   leider müssen wir hier eine kleine Toleranz addieren,
        %       da Matlab 2.0 >= oder <= 2 nicht als boolsches Wahr interpretiert,
        %       um jedoch die Verfälschung möglichst klein zu halten,
        %       müssen wir auch den addierten Term klein wählen)
        if norm(d) > tol && all(([G;U] * (x+d)) <= ([b r]+1.0e-10)')
            x = x + d;
            continue;
        end
        %Case 4
        %   leider tritt hier der gleiche Fall wie im Case 3 auf.
        if norm(d) > tol && all(([G;U] * (x+d)+1.0e-10) >= ([b r])')
            u_i = U;
            u_i(W,:) = [];
            r_i = r;
            r_i(W) = [];
            t = [];
            for i = 1:size(u_i,1)
                t(i) = (r_i(i)-u_i(i,:)*x)/(u_i(i,:)*d);
            end
            [t_min,min_index_t] = min(t);
            x = x+t_min*d;
            W(end+1) =  min_index_t;
            continue;
        end      
    end 
end
end

