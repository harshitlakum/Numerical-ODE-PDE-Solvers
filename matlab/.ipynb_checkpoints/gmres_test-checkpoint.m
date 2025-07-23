% gmres_test.m
% Test restarted GMRES on two unsymmetric matrices, with and without ILU

function gmres_test
    %% Parameters
    restarts = [10, 20, 40, 80];
    tol      = 1e-6;
    
    %% Build test problems
    problems = {
      'Convection-Diffusion (50x50)', convection_diffusion_2d(50,1.0,50.0);
      'Random DD Sparse (5000x5000)', build_random_dd(5000,1e-3)
    };
    
    %% Loop over problems
    for k = 1:size(problems,1)
        name = problems{k,1};
        A    = problems{k,2};
        b    = ones(size(A,1),1);
        
        fprintf('\n=== %s ===\n', name)
        fprintf(' Precond    m    Iterations\n')
        fprintf(' -------------------------\n')
        
        % Attempt ILU
        try
            setup.type    = 'ilutp';
            setup.droptol = 1e-4;
            [L,U] = ilu(A, setup);
            precs = {'None',[]; 'ILU',{L,U}};
        catch
            precs = {'None',[]};
        end
        
        % For each preconditioner and restart
        for p = 1:size(precs,1)
            preName = precs{p,1};
            if isempty(precs{p,2})
                M1 = []; M2 = [];
            else
                M1 = precs{p,2}{1};
                M2 = precs{p,2}{2};
            end
            
            for m = restarts
                [~, flag, relres, iterInfo] = ...
                  gmres(A, b, m, tol, 1, M1, M2);
                totalIters = iterInfo(1)*m + iterInfo(2);
                fprintf(' %-8s %4d %12d\n', preName, m, totalIters)
            end
        end
    end
end

function A = convection_diffusion_2d(n, diffusion, convection)
    N = n*n;
    e = ones(N,1);
    main      = 4*diffusion * e;
    ew_diff   = -diffusion * ones(N-1,1);
    ns        = -diffusion * ones(N-n,1);
    bias      = convection/(2*(n-1));
    ew_conv_p = -bias * ones(N-1,1);
    ew_conv_m =  bias * ones(N-1,1);
    ew_plus   = ew_diff + ew_conv_p;
    ew_minus  = ew_diff + ew_conv_m;
    diags = [main, ew_minus, ew_plus, ns, ns];
    offs  = [  0,      -1,      +1, -n, +n];
    A = spdiags(diags, offs, N, N);
    for i = 1:(n-1)
        idx = i*n;
        A(idx,   idx-1) = 0;
        A(idx-1, idx  ) = 0;
    end
    A = sparse(A);
end

function A = build_random_dd(n, density)
    R1 = sprand(n, n, density, 1);
    R2 = sprand(n, n, density, 2);
    A  = R1 - R2 + spdiags(10*ones(n,1), 0, n, n);
end