clc;
clear all;

a1 = 2;
a2 = 3;
b1 = 3.9;
b2 = 1;
h = 1;

n = 12;
A = zeros(n,n);

%C
for i = 1:n
    A(i,i) = 2*a1+2*a2;
end

%E
for i = 1:n-1
    if i == 4 || i ==8
        A(i,i+1) = 0;
    else
        A(i,i+1) = (b1*h)/2 - a1;
    end
end

%W    
for i = 1:n-1
    if i == 4 || i ==8
        A(i,i+1) = 0;
    else
        A(i+1,i) = (-a1)-(b1*h)/2;
    end
end

%N
for i = 1:n-4
    A(i,i+4) = (b2*h)/2 - a2;
end

%S
for i = 1:n-4
    A(i+4,i) = (-b2*h)/2 - a2;
end

%set linear systems Ax=b
b = zeros(n,1);
for i = 1 : n/2
    b(i,1) = 1;
end

%gmres

tol = 1e-12;
maxit = 12;
[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit);

gm = cond(A) %condtion number of A when doing gmres solving the Ax=b 

%ilu preconditioned gmres
B = sparse(A);
[L,U] = ilu(B,struct('type','ilutp','droptol',1e-6));
[x1,fl1,rr1,it1,rv1] = gmres(B,b,[],tol,maxit,L,U);
lu=full(L)*full(U);
precond_gm = cond((lu)^(-1)*B) %condition number after add ilu preconditioner

hold on
%plot gmres
semilogy(0:maxit,rv0/norm(b),'-o');
dim = [.35 .22 .3 .3];
str = {'condition number: 6.0365'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%plot preconditioned gmres
semilogy(0:it1(2),rv1/norm(b),'-*');
dim = [.19 .007 .3 .18];
str = {'condition number: 1.0000'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

xlabel('Iteration number');
ylabel('Relative residual');
legend('GMRES (without restart and preconditioned)','iLU-Preconditioned GMRES');
