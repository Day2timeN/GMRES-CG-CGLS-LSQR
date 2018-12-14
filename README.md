# GMRES-CG-CGLS-LSQR
The program breaks into three parts. First part is to use convection-diffusion test matrices to test the effect of ILU preconditioned on GMRES convergence. 

Second part is to test the numerical differences on CG, CGLS and LSQR by solving normal equation AT Ax = AT b. Also, choosing A by selecting α on L−αI, where L is the finite-difference Laplacian matrix and I is identity. 

Last part is to find a good preconditioner and solve the least square problem by LSQR.
