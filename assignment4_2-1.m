clc;
clear all;

%construct A
A = zeros(4,4);
for i = 1:4
    A(i,i) = 4;
end
for i = 1:3
    A(i,i+1) = -1;
end
for i = 1:3
    A(i+1,i) = -1;
end

n = 3; %12 by 12 matrix
Ar = repmat(A,1,n);
Ac = mat2cell(Ar,size(Ar,1),repmat(size(A,1),1,n));
diag_A = blkdiag(Ac{:});

%construct -1 bond of A
A1 = zeros(4*n,4*n);
for i = 1 : 4*n-4
    A1(i,i+4) = -1;
end
for i = 1 : 4*n-4
    A1(i+4,i) = -1;
end

A = A1 + diag_A; %result A

%construct b
b = zeros(4*n,1);
for i = 1 : 2*n
    b(i) = 1;
end

%find A 
a = linspace(-20,20,1000);
cond_A = zeros(length(a),1);
cond_At = zeros(length(a),1);
for i = 1:length(a)
    A_test = A - a(i)*eye(4*n); %L - alpha*I
    cond_A(i) = cond(A_test);
    cond_At(i) = cond((A_test)'*A_test);
end

hold on
title('Analysis of Condition Number for 12 by 12 Matrix');
plot(a,cond_A,'-*');
plot(a,cond_At,'-o');
ylim([0 1000]);
xlabel('alpha');
ylabel('condition number')
legend('Condition number of A','Condition number of A^TA');





