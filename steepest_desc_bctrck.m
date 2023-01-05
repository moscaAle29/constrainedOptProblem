function [xk,fk,gradf_norm, k, xseq,btseq,deltaxk_norm]=constr_steepest_desc_bctrck...
    (x0,f,gradf,kmax,tolgrad,c1,rho,btmax,gamma, tolx,pi_X)

    xk=x0;
    k=0;
    xseq=zeros(size(x0));
    xseq(:,1)=x0;
    reached=0;
    btseq=zeros(kmax);
    for k=1:kmax
        xkm1=xk;
        alphak=1;
        count=0;
        p=-gradf(xk);
        xbar=xk+gamma*p;
        xhatk=pi_X(xbar);
        pik=(xhatk-xk);
        xk=xk+alpha*pik;
        while count<=btmax && f(xk+alphak*(pik))>armijo(xk,c1,alphak,f,gradf)
            alphak=alphak*rho;
            count=count+1;
        end
        btseq(k)=count;
        p=-gradf(xk);
        xseq=[xseq xk];
        gradf_norm=norm(p);
        deltaxk_norm=norm(xkm1-xk);
        if norm(p)< tolgrad || deltaxk_norm<=tolx
           reached=1;
           fprintf('reached the tollerance %f\n',tolgrad)
           break
        end
    end
    fk=f(xk);
    x_seq=xseq(:,k);
    if reached==0 
        fprintf('couldnt reach the tollerance %d at step %d\n',tolgrad,kmax)
    end
end

function value=armijo(xk,c1,ak,f,gradf)
    value=f(xk)+c1*ak*norm(gradf(xk));
end