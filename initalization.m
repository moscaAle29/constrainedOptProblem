clear all
close all
clc

%define the variables for the problem
d=5;
n=10^d;
x0=100*ones(n,1);
alpha0=1;
kmax=1000;
tolgrad=10^-3;
c1=10^-4;
rho=0.8;
btmax=50;
gamma=0.1;
tolx=10^-12;

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

x_prime = zeros(n,1);
for i = (n/2 + 1) : n
    x_prime(i) = 1;
end
