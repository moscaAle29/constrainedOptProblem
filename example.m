initalization;

%define handles
Pi_X = @(x) box_projection(x, X3);
f = @(x) f(x);
fin_gradf = @(x) findiff_grad(f,x,10^-4,'fw');

%adjust tolgrad
tolgrad = tolgrad + norm(x_prime);

[xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
    prof_constr_steepest_desc_bcktrck(x0, f, fin_gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);









