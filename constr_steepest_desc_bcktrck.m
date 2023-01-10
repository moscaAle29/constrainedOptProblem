function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X)

xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = Pi_X(x0);
fk = f(xk);
gradfk = gradf(xk);

k = 0;
gradfk_norm = norm(gradfk);
deltaxk_norm = tolx + 1;

farmijo = @(fk, alpha, gradfk, pk) fk + c1 * alpha * gradfk' * pk;

while k < kmax && gradfk_norm >= tolgrad && deltaxk_norm >= tolx
    
    pk = -gradf(xk);
    
    xbark = xk + gamma * pk;
    xhatk = Pi_X(xbark);    
    
    alpha = 1;
    
    pik = xhatk - xk;
    xnew = xk + alpha * pik;
    fnew = f(xnew);
    
    bt = 0;
    
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pik)
        
        alpha = rho * alpha;

        xnew = xk + alpha * pik;
        fnew = f(xnew);

        bt = bt + 1;
        
    end
    
    deltaxk_norm = norm(xnew - xk);
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    
    k = k + 1;
    
    xseq(:, k) = xk;
    btseq(k) = bt;
end

xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end