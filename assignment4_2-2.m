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

s = 4*n;
r = zeros(s,1);
u_CG = zeros(s,1);
u_CG(1:s/2,1) = 1; %initial approximation
p = zeros(s,s); %direction
a = zeros(s,1); %step size
beta = zeros(s,1); %beta
temp_AP = zeros(s,s); %store A*P
temp_rr = zeros(s,1); %store (r,r)

tol2 = 10^(-8);
relerr2 = inf;
niter2 = 1;

%solve normal equation ATAx = ATb by CG method
A_L = A -5*eye(4*n); %A = L - alpha*I
b_t = A_L'*b; %AT*b
A_t = A_L'*A_L; %AT*A
cond(A_t)
eig(A_t)
u_exact = A_t^(-1)*b_t;

r(:,1) = A_L'*(b - A_L*u_CG(:,1)); %r0
p(:,1) = r(:,1); %p0

while relerr2 > tol2
    % store it and use it after in order to do only once matrix-vector product
    temp_AP(:,niter2) = A_t*p(:,niter2); 
    temp_rr(niter2,1) = dot(r(:,niter2),r(:,niter2)); % do same thing above in inner product
    a(niter2,1) = (dot(r(:,niter2),r(:,niter2)))/(dot((temp_AP(:,niter2)),p(:,niter2))); %new alpha
    u_CG(:,niter2+1) = u_CG(:,niter2) + a(niter2,1)*p(:,niter2); %new u
    
    %updates parameters
    r(:,niter2+1) = r(:,niter2) - a(niter2,1)*temp_AP(:,niter2); %update r
    beta(niter2,1) = (dot(r(:,niter2+1),r(:,niter2+1)))/temp_rr(niter2,1); %update beta
    p(:,niter2+1) = r(:,niter2+1) + beta(niter2,1)*p(:,niter2); %update direction p
    
    relerr2 = norm(b_t - A_t*u_CG(:,niter2));
    niter2 = niter2 + 1;
end

%residual
residual_CG = zeros(niter2,1);
for i = 1:niter2
    residual_CG(i) = norm(b_t - A_t*u_CG(:,i));
end

%CGNR or CGLS
r1 = zeros(s,1);
z = zeros(s,1);
u_CGLS = zeros(s,1);
u_CGLS(1:s/2,1) = 1; %initial approximation
r1(:,1) = b - A_L*u_CGLS(:,1); %r0
z(:,1) = A_L'*r1(:,1); %z0
p1 = zeros(s,s);
p1(:,1) = z(:,1); %p0
q = zeros(s,s); %q0
a1 = zeros(s,1); %alpha
beta1 = zeros(s,1); %beta

tol1 = 10^(-8);
relerr1 = inf;
niter1 = 1;

while relerr1 > tol1
    q(:,niter1) = A_L*p1(:,niter1); %q = A*p0
    %alpha_k = s_k^T*s_k/q_k^T*q_k
    a1(niter1,1) = (z(:,niter1)'*z(:,niter1))/(q(:,niter1)'*q(:,niter1)); 
    u_CGLS(:,niter1+1) = u_CGLS(:,niter1) + a1(niter1,1)*p1(:,niter1); %x_k+1 = x_k + alpha_k*p_k
    r1(:,niter1+1) = r1(:,niter1) - a1(niter1,1)*q(:,niter1); %r_k+1 = r_k - alpha_k*q_k
    z(:,niter1+1) = A_L'*r1(:,niter1+1); %z_k+1 = A^T*r_k+1
    beta1(niter1,1) = (z(:,niter1+1)'*z(:,niter1+1))/(z(:,niter1)'*z(:,niter1));%beta_k = s_k+1^T*s_k+1/s_k^T*s_k
    p1(:,niter1+1) = z(:,niter1+1) + beta1(niter1,1)*p1(:,niter1);%p_k+1 = s_k+1 + beta_k*p_k
    
    relerr1 = norm(b-A_L*u_CGLS(:,niter1));
    niter1 = niter1 + 1;
end

%residual
residual_CGLS = zeros(niter1,1);
for i = 1:niter1
    residual_CGLS(i) = norm(b-A_L*u_CGLS(:,i));
end

%LSQR
tol = 10^(-8);
maxit = 17;
M1 = [];
M2 = [];
x0 = zeros(s,1);
x0(1:size(A,1)/2,1) = 1;
[x,flag,relres,iter,resvec] = lsqr(A_t,b_t,tol,maxit,M1,M2,x0);


hold on
plot(residual_CG,'-*');
plot(residual_CGLS,'-o');
plot(resvec,'-x');

title('Analysis of convergence when A^TA is not well-conditioned, condition number ~= 391');
dim = [.5 .22 .3 .3];
str = {'For alpha = 3'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlabel('number of iterates');
ylabel('residual')
legend('CG','CGLS','LSQR');

