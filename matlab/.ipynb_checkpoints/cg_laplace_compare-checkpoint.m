% cg_laplace_compare.m
% CG for 2D Laplace with/without incomplete Cholesky; convergence plots

n = 100;                   % n^2 = system size
N = n*n;
e = ones(n,1);
K1d = spdiags([-e 2*e -e], -1:1, n, n);
I = speye(n);
A = kron(I, K1d) + kron(K1d, I);

x_true = rand(N,1);        % Random exact solution
b = A*x_true;
x0 = zeros(N,1);

resvec_noprec = [];
[x1,flag1,relres1,iter1,resvec_noprec] = pcg(A, b, 1e-8, 1000, [], [], x0);

try
    L = ichol(A, struct('droptol',1e-4,'type','ict'));
    resvec_prec = [];
    [x2,flag2,relres2,iter2,resvec_prec] = pcg(A, b, 1e-8, 1000, L, L', x0);
catch err
    warning('Incomplete Cholesky failed: %s', err.message);
    resvec_prec = [];
end

figure; hold on;
semilogy(0:length(resvec_noprec)-1, resvec_noprec, 'b-', 'LineWidth',1.5)
if ~isempty(resvec_prec)
    semilogy(0:length(resvec_prec)-1, resvec_prec, 'r-', 'LineWidth',1.5)
    legend('CG (no preconditioner)','CG (incomplete Cholesky)','Location','northeast')
else
    legend('CG (no preconditioner)','Location','northeast')
end
xlabel('Iteration')
ylabel('Residual Norm')
title(sprintf('CG Convergence: 2D Laplacian, N=%d',N))
grid on; box on;