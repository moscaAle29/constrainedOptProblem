function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx,...
    findiff_enable, k_findiff,method,feasible_set)
%{
=================================================================
input:
- xo -> starting point for the method
- f -> function handle for the function to optimize
- gradf -> function handle for the exact gradient
- kmax -> max number of iterations for the method
- tolgrad -> minimum increase in th gradient to be considered 0 and end the
            iterations
- c1 -> constant for armijo
- rho -> coefficient for backtrack
- btmax -> max number of iteration for the backtrack
- gamma -> discount factor if we are out of the feasible set
- tolx -> minimum increase in the value x to be considered 0
- findiff_enable -> enable the finite difference computation of the gradient

output:
- xk -> x value at the last iteration
- fk -> f value at the last iteration
- gradfk_norm -> norm of the gradient at the last iteration
- deltaxk_norm -> difference between xk and the new direction
- k -> last iteration of the method
- xseq -> sequence of points the method goes through
- btseq -> sequence of points in the backtracking method
=================================================================
%}
geth=@(xhat,k) norm(xhat)*10^-k;

xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = box_projection(x0,feasible_set);
fk = f(xk);
if findiff_enable
    gradfk = findiff_grad(f,xk,get_h(xk,k_findiff),method);
else
    gradfk = gradf(xk);
end
k = 0;
gradfk_norm = norm(gradfk);
deltaxk_norm = tolx + 1;

farmijo = @(fk, alpha, gradfk, pk) fk + c1 * alpha * gradfk' * pk;

while k < kmax && gradfk_norm >= tolgrad && deltaxk_norm >= tolx
    
    if findiff_enable
        pk = -findiff_grad(f,xk,geth(xk,k_findiff),metod);
    else
        pk = -gradf(xk);
    end
    
    xbark = xk + gamma * pk;
    xhatk = box_projection(xbark,feasible_set);    
    
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
    if findiff_enable
        gradfk = -findiff_grad(f,xk,get_h(xk,k_findiff),method);
    else
        gradfk = -gradf(xk);
    end
    gradfk_norm = norm(gradfk);
    
    k = k + 1;
    
    xseq(:, k) = xk;
    btseq(k) = bt;
end

xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end
