% =====================================
% Constrained Optimization Porject
% 
% Alessandro Mosca -> s309595
% Giacomo Qualcosa -> s
% Thriet Min Ngo -> s309062
% 
% projected Gradient Method
% ====================================
%define the variables for the problem
d=3;
n=10^3;
%note for Triet and Giacomo, i think we should run it 3 times, with d=3, d=4, d=5 
%start by defining the feasible set
X1=zeros(n,2);
X2=zeros(n,2);
X3=zeros(n,2);

%define the first feasible set
for i=1:n
    X1(i,1)=1;
    X1(i,2)=5.12;
end

%define the second feasible set
X2(1,1)=-5.12;
X2(1,2)=5.12;
for i=2:n
    X2(i,1)=1;
    X2(i,2)=5.12;
end

%define the third feasible set
%note for Giacomo and Triet: doing this is computationally more efficient
%because the processor doesn't have to compute the conditional jump of the
%if, if you want more clarity we can impmlement it in a different way
for i=1:n/2
    X3(i,1)=-5.12;
    X3(i,2)=5.12;
end

for i=(n/2+1):n
    X3(i,1)=1;
    X3(i,2)=5.12;
end


%solve the unconstrained problem with exact gradient 
[xk_u_e,fk_u_e,gradf_norm_u_e, k_u_e, xseq_u_e,btseq_u_e]=steepest_desc_bctrck(x0,...
    f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax,false);

%solve the unconstrained problem with finite differences of the gradient
for findiff_k=1:2:12
    [xk_u_f,fk_u_f,gradf_norm_u_f, k_u_f, xseq_u_f,btseq_u_f]=...
        constr_steepest_desc_bcktrck(x0,f,gradf,kmax,tolgrad, c1,rho,btmax,...
        gamma,tolx,findiff_k,true);
end
%solve the constrained problem with exact gradient
[xk_c_e, fk_c_e, gradfk_norm_c_e, deltaxk_norm_c_e, k_c_e, xseq_c_e, btseq_c_e] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf,kmax, tolgrad, c1, rho,...
    btmax, gamma, tolx,k_findiff,false);

%solve the constrained problem with finite differences
for findiff_k=1:2:12
    [xk_c_f, fk_c_f, gradfk_norm_c_f, deltaxk_norm_c_f, k_c_f, xseq_c_f, btseq_c_f] = ...
        constr_steepest_desc_bcktrck(x0, f, gradf,kmax, tolgrad, c1, rho,...
        btmax, gamma, tolx,k_findiff,false);
end
% define a function handle for the function
function f_value = f(x)
    sum=0;
    len=length(x);
    for i= 1:len
        sum=sum+i*x(i)^2;
    end
    f_value=sum;
end

%define a function handle fot the exact gradient
%note for Triet and Giacomo am i dumb for the computation of the gradient?
function gradf_value= gradf(x)
    len_x=length(x);
    gradf_value=zeros(len_x);
    for i=1:len_x
        gradf_value(i)=2*i*x(i);
    end
end