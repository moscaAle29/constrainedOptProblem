function [xk, fk, gradfk_norm, k, xseq, btseq] = steepest_desc_bcktrck(x0, ...
    f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax,findiff_enable)

xk = x0;
fk = f(xk);
gradfk = extractdata(gradf(xk));
gradfk_norm = norm(extractdata(gradf(xk)));
k = 0;
xseq = zeros(length(x0),kmax); %meglio così per evitare di occupare troppa 
                               %memoria (se lo inizializzassi come un 
                               %vettore vuoto ogni volta dovrei
                               %sovrascrivere
btseq = zeros(1, kmax);

Armijo = @(f, alpha, gradient, p) f + c1 * alpha * gradient' * p;

    while k < kmax && gradfk_norm >= tolgrad
        
        alphak = alpha0;
        pk = -gradf(xk);
        xnew = xk + alphak * pk;
        fnew = f(xnew);
        bt = 0;

        while bt <= btmax && fnew > Armijo(fk, alphak, gradfk, pk)

            alphak = rho * alphak;
            xnew = xk + alphak * pk;
            fnew = f(xnew);
            bt = bt + 1;

        end

        xk = xnew;
        fk = fnew;
        gradfk = extractdata(gradf(xk));
        gradfk_norm = norm(gradfk);

        k = k+1;
        xseq(:,k) = xk;
        btseq(k) = bt;
        
    end

    fk = f(xk);
    xseq = [x0, xseq(:, 1:k)];
    btseq = btseq(1:k);

end