B = full(A);
left = B'*b;
C = B'*B;
C1 = sparse(C);

%use incomplete cholesky with shifted to find preconditioner for A^TA
%since "pivot breakdown", we need to increase the diagonal dominance of A^T
%A by diagonal shifts
alpha = 1.3;
L1 = ichol(C1, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
[x,flag,relres,iter,resvec] = lsqr(B,b,1e-6,10000,L1,[]);

%LSQR without preconditioner
[x1,flag1,relres1,iter1,resvec1] = lsqr(B,b,1e-6,10000,[],[]);

%QR factorization by householder algorithm
[Q,R]=qr(B,0); %lots of fill-in
[x2,flag2,relres2,iter2,resvec2] = lsqr(B,b,1e-6,10000,R,[]);

%used modified incomplete cholesky with no fill-in
L2 = ichol(C1, struct('type','nofill','michol','on','diagcomp',1.2));
[x3,flag3,relres3,iter3,resvec3] = lsqr(B,b,1e-6,10000,L2,[]);

%QR factorization obtained from modified Gram-Schmidt 
G = B;
[m,n] = size(G);
r = zeros(n,n);
q = zeros(m,n);
for i = 1 : n
    r(i,i) = norm(G(:,i));
    q(:,i) = G(:,i)/r(i,i);
    for j = i+1 : n
        if r(i,j) ~= 0
            r(i,j) = q(:,i)'*G(:,j);
            G(:,j) = G(:,j) - r(i,j)*q(:,i);
        end
    end
end
[x4,flag4,relres4,iter4,resvec4] = lsqr(B,b,1e-6,10000,r,[]);

%QR factorization obtained from incomplete modified Gram-Schmidt
%sparsity pattern
S = B;
r1 = zeros(n,n);
q1 = zeros(m,n);
for i = 1:n
    r1(i,i) = norm(S(:,i));
    q1(:,i) = S(:,i)/r1(i,i);
    for j = i+1:n
        if B(i,j)~= 0
            r1(i,j) = q1(:,i)'*S(:,j);
        else
            r1(i,j) = 0;
        end
        S(:,j) = S(:,j) - r1(i,j)*q1(:,i);
    end
end
[x5,flag5,relres5,iter5,resvec5] = lsqr(B,b,1e-6,10000,r1,[]);

hold on
xlim([0 20])
xlabel('number of iteration')
ylabel('residual')
plot(resvec,'-*')
plot(resvec1,'-o')
plot(resvec2,'-x')
plot(resvec3,'-+')
plot(resvec4,'--')
plot(resvec5,'-v')
legend('incomplete-cholesky preconditioner','no preconditioner','QR-factorization(Household) preconditioner','no-fill-in modified-incomplete-cholesky preconditioner','QR-factorization(modified Gram-Schmidt) preconditioner','QR-factorization(incomplete modified Gram-Schmidt) preconditioner')








