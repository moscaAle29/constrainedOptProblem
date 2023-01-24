% =====================================
% Constrained Optimization Porject
%
% Bastiani Giacomo -> s303217
% Mosca Alessandro -> s309595 
% Ngo Minh Triet-> s309062
% 
%projected Gradient Method
% ====================================
clear all
close all
%define the variables for the problem
d=5;
n=10^d;
x0=100*ones(n,1);
alpha0=1;
kmax=1000;
tolgrad=10^-(8-d);
c1=10^-4;
rho=0.8;
btmax=50;
gamma=0.1;
tolx=10^-(8-d);

%start by defining the feasible set
X1=zeros(n,2);
X2=zeros(n,2);
X3=zeros(n,2);

% define a function handle for the function
f=@(x) sum((1:length(x))'.*(x.^2));
%define a funcion handle for the gradient of the function
gradf=@(x) 2*(1:length(x))'.*x ;

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
tic
[xk,fk,gradf_norm, k, xseq,btseq]=...
    steepest_desc_bcktrck(x0,f, gradf, alpha0, kmax,...
    tolgrad, c1, rho, btmax,false,0,'fw');
unconstrained_time_e=toc*ones(6,1);
unconstrained_error_e=norm(xk)*ones(6,1);

%solve the unconstrained problem with finite differences of the gradient
%using forward method
unconstrained_time_findiff_fw=zeros(6,1);
unconstrained_error_findiff_fw=zeros(6,1);
i=1;
for findiff_k=2:2:12
    tic;
    [xk, fk, grad_norm, k, xseq,btseq]=...
        steepest_desc_bcktrck(x0,f,gradf,alpha0 ,kmax, tolgrad,...
        c1, rho, btmax,true,findiff_k,'fw'); %#ok<*ASGLU> 
    unconstrained_time_findiff_fw(i)=toc;
    unconstrained_error_findiff_fw(i)=norm(xk);
    i=i+1;
end
%solve the unconstrained problem with finite differences of the gradient
%using centered method
unconstrained_time_findiff_cd=zeros(6,1);
unconstrained_error_findiff_cd=zeros(6,1);
i=1;
for findiff_k=2:2:12
    tic;
    [xk,fk,gradf_norm, k, xseq,btseq]=...
        steepest_desc_bcktrck(x0,f,gradf,alpha0 ,kmax, tolgrad,...
        c1, rho, btmax,true,findiff_k,'c');
    unconstrained_time_findiff_cd(i)=toc;
    unconstrained_error_findiff_cd(i)=norm(xk);
    i=i+1;
end
%plot the time it takes us to execute the function
x=linspace(1,6,6);
%plot the time
figure(1)
plot(x,unconstrained_time_e,x,unconstrained_time_findiff_fw,x,unconstrained_time_findiff_cd)
title('time used to solve unconstrained problem')
xlabel("value of k")
ylabel("seconds")
legend({'exact gradient time','forward finite differences time','centered finite differences time'},'location','southeast')
%plot the error
figure(2)
plot(x,unconstrained_error_e,x,unconstrained_error_findiff_fw,x,unconstrained_error_findiff_cd)
title('norm of the error in unconstrained problem')
xlabel("value of k")
ylabel("error")
legend({'exact gradient error','forward finite differences gradient error','centered finite differences gradient error'},'location','southeast')

%set the feasible set
feasible_set=X2;
Pi_X = @(x) box_projection(x, X3);
%solve the constrained problem with exact gradient 
tic
[xk,fk,gradf_norm,deltaxk_norm, k, xseq,btseq]=...
    constr_steepest_desc_bcktrck(x0,f, gradf, kmax,...
    tolgrad, c1, rho, btmax, gamma, tolx, feasible_set, ...
    false, 0,'fw');
constrained_time_e=toc*ones(6,1);
constrained_error_e=norm(xk-box_projection(zeros(n,1),feasible_set))*ones(6,1);

%solve the constrained problem with finite differences of the gradient
%using forward method
constrained_time_findiff_fw=zeros(6,1);
constrained_error_findiff_fw=zeros(6,1);
i=1;
for findiff_k=2:2:12
    tic;
    [xk,fk,gradf_norm,deltaxk_norm, k, xseq,btseq]=...
    constr_steepest_desc_bcktrck(x0,f, gradf, kmax,...
    tolgrad, c1, rho, btmax, gamma, tolx, feasible_set,...
    true, findiff_k,'fw');
    constrained_time_findiff_fw(i)=toc;
    constrained_error_findiff_fw(i)=norm(xk-box_projection(zeros(n,1),feasible_set));
    i=i+1;
end
%solve the unconstrained problem with finite differences of the gradient
%using centered method
constrained_time_findiff_cd=zeros(6,1);
constrained_error_findiff_cd=zeros(6,1);
i=1;
for findiff_k=2:2:12
    tic;
    [xk,fk,gradf_norm,deltaxk_norm, k, xseq,btseq]=...
    constr_steepest_desc_bcktrck(x0,f, gradf, kmax,...
    tolgrad, c1, rho, btmax, gamma, tolx, feasible_set,...
    true, findiff_k,'fw');
    constrained_time_findiff_cd(i)=toc;
    constrained_error_findiff_cd(i)=norm(xk-box_projection(zeros(n,1),feasible_set));
    i=i+1;
end
%plot the time it takes us to execute the function
x=linspace(1,6,6);
%plot the time
figure(3)
plot(x,constrained_time_e,x,constrained_time_findiff_fw,x,constrained_time_findiff_cd)
title('time used to solve constrained problem')
xlabel("value of k")
ylabel("seconds")
legend({'exact gradient time','forward finite differences time','centered finite differences time'},'location','southeast')
%plot the error
figure(4)
plot(x,constrained_error_e,x,constrained_error_findiff_fw,x,constrained_time_findiff_cd)
title('norm of the error in constrained problem')
xlabel("value of k")
ylabel("error")
legend({'exact gradient error','forward finite differences gradient error','centered finite differences gradient error'},'location','southeast')
