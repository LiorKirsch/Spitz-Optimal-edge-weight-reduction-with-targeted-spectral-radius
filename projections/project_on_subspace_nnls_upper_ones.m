function [edge_weights,vertices_weights] = project_on_subspace_nnls_upper_ones(b,A,AtA,x0)
% project_on_subspace_nnls solves nonnegative least squares with each
% variable is also bounded above by 1
% It uses the design matrix D (and the precomputed DtD) to solve nnls using
% lbfgs.
% For performance we use the LBFGSB package
% 
% [edge_weights,vertices_weights] = project_on_subspace_nnls(x_values,D,DtD);
%

    [num_edges,N] = size(A);
    
    fcn     = @(x) norm( A*x - b)^2;
    % here are two equivalent ways to make the gradient. grad2 is sometimes faster
%     grad1    = @(x) 2*A'*(A*x-b);

    Ab = A'*b;
    grad2    = @(x) 2*( AtA*x - Ab );

    grad    = grad2;

    l  = zeros(N,1);    % lower bound
    u  = ones(N,1);      % ones upper bound 

    fun     = @(x)fminunc_wrapper( x, fcn, grad); 
    % Request very high accuracy for this test:
    % m is the number of lbfgs vectors
    % factr - Tolerance
    % pgtol - Another tolerance setting, relating to norm(gradient,Inf)
    % printEvery - print every X iterations.

    opts    = struct( 'factr', 1e2, 'pgtol', 1e-10, 'm', 10);
    opts.printEvery  = inf;
    
    if exist('x0','var')
        opts.x0 = x0;
    end
    
    if N > 10000
        opts.m  = 50;
    end

    [vertices_weights, ~, info] = lbfgsb(fun, l, u, opts );
   
    edge_weights = A*vertices_weights;
end