% =====================================
% Constrained Optimization Porject
% 
% Alessandro Mosca -> s309595
% Giacomo Qualcosa -> s
% Thriet Min Ngo -> s
% 
% projected Gradient Method
% ====================================
%define the variables for the problem
d=3;
n=10^3;

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

% define a function handle for the function
function y = f(x)
    sum=0;
    len=length(x);
    for i= 1:len
        sum=sum+i*x(i);
    end
    y=sum;
end

%solve the unconstrained problem 

%solve the constrained problem

%{

for i = 1:2:12
    h=norm(xhat)*10^-i;
end
%}
