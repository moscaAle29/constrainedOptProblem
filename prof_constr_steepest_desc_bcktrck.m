function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
    prof_constr_steepest_desc_bcktrck(x0, f, gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X)
%
% function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
%     constr_steepest_desc_bcktrck(x0, f, gradf, alpha0, ...
%     kmax, tollgrad, c1, rho, btmax, gamma, tolx, Pi_X)
%
% Projected gradient method (steepest descent) for constrained optimization.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
% gamma = the initial factor that multiplies the descent direction at each
% iteration;
% tolx = value used as stopping criterion w.r.t. the norm of the
% steps;
% Pi_X = projection function
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% deltaxk_norm = length of the last step of the sequence
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = Pi_X(x0); % Project the starting point if outside the constraints
fk = f(xk);
gradfk = gradf(xk);

k = 0;
gradfk_norm = norm(xk);
deltaxk_norm = tolx + 1;

while k < kmax && gradfk_norm >= tolgrad && deltaxk_norm >= tolx
    disp(k)
    % Compute the descent direction
    pk = -gradfk;
    
    xbark = xk + gamma * pk;
    xhatk = Pi_X(xbark);    
    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    pik = xhatk - xk;
    xnew = xk + alpha * pik;
    
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo (w.r.t. pik) condition not satisfied
    
    
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pik)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pik;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm, deltaxk_norm
    deltaxk_norm = norm(xnew - xk);
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(xk);
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end