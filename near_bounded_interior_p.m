function [W,exitflag] = near_bounded_interior_p(A, beta, options)
%     This function computes the nearest bounded matrix to the matrix A.
%     Each element in the output matrix W is bounded above by the value of A
%     and bounded bellow by 0.
%     The largest eigenvalue of W is bounded by beta.
%     If the matrix A is symmetric we only use half the variables. 
%
%     Input: 
%           A  - sparse adjacancy matrix.
%           beta - the bound for the first eigen value (the spectral norm)
%
%           options.Objective - choose the objective L2 (default) or L1.
%           options.TolFun    - is a convergence tolerance, which defaults to 16*EPS. 
%           options.MaxIter - the maximum number of iterations (default 100, but may
%               need to be increased, depends on the number of eigen values that
%               we need to cap).
%           options.Display   - a boolean flag which contronls the print
%               after infinit number of steps, we use a small dist tp find a point 
%               inside the convex group of matrices with bounded spectrum
%               (default 1e-10)
%           options.UseRadius - if A is not symmetric we can use the svd (default) or
%               the schur factorization to get the bounds on the radius
%               true - specifices using schur, false (default) uses svd
%           
%     Example: 
%           dim = 100;       num_nnz = dim^2 /4;       beta = 10;
% 
%           row_indices = randi(dim,num_nnz,1);
%           column_indices = randi(dim,num_nnz,1);
%           values = rand(num_nnz,1);
%           A = sparse(row_indices, column_indices,values, dim,dim);
%           A = A - diag(diag(A));
%           A = A + A';
% 
%           W = near_bounded_interior_p(A, beta, []);
%
if ~exist('options','var') options = struct; end
if isfield(options,'use_my_svd') 
    use_my_svd=options.use_my_svd;
else
    use_my_svd=false;
end
if isfield(options,'UseRadius') 
    use_radius=options.UseRadius;
else
    use_radius=false;
end
if isfield(options,'Objective') 
    objective=options.Objective;
else
    objective='L2';
end

    
        
    % ==== Perform some checks and statistics ====
    [dim, dim2] = size(A);
    assert(dim==dim2,'A should be a squared matrix');
    
    is_sparse = issparse(A);
    issym = issymmetric(A);
%     issym = true ; 
    
    if issym
        fprintf('using symmetric properties\n');
    end
    opts.issym = issym;

    top2eigs = eigs(A,2,'lm', opts);
    start_2norm = normest(A);
    fprintf('Spectral radius of the original problem is %g\n', abs(top2eigs(1)));
    fprintf('Spectral  norm  of the original problem is %g\n', start_2norm );
    fprintf('Spectral gap is  %g\n', ( abs(top2eigs(1)) - abs(top2eigs(2)) ));

%     if is_sparse
        if issym
            [row_indices, column_indices, values] = extract_upper_half_inds(A);
        else
            [row_indices, column_indices, values] = find(A);
        end
%     else
%         values = A(:);
%         [column_indices,row_indices] = meshgrid(1:dim,1:dim);
%         row_indices = row_indices(:);
%         column_indices = column_indices(:);
%     end
   
    % ===   using the fact that ||A||_2 <= ||A||_F or that R(A) <= ||A||_F ===
    x0 = values * min(1, beta/(norm(values) + eps));
    
    % ==== Some optimization parameters ====
    defoptions = optimset('Algorithm','interior-point');
    defoptions.Display = 'iter' ;%'final'; %'iter';
    defoptions.MaxFunEvals = 5000;
    defoptions.GradObj = 'on';
    defoptions.GradConstr = 'on';
    defoptions.Hessian = {'lbfgs',6};

%     defoptions.AlwaysHonorConstraints = 'none';
%     defoptions.SubproblemAlgorithm = 'cg';
%     defoptions.DerivativeCheck = 'on';
%     defoptions.ScaleProblem = 'obj-and-constr';
%     defoptions.TolFun = 1e-14;
    options = optimset(defoptions, options);

    switch objective
        case 'L1'
            obective_func = @(x) func_objective_L1(x,values);
        case 'L2'
            obective_func = @(x) func_objective(x,values);
        otherwise
            error('unkown objective - %s, should be ''L1'' or ''L1'' ', objective);
    end
    
    lb = zeros(size(values));
    ub = values; % before removing this be sure to check the L1 objective
    % ==== Solve the optimization problem ====
    [x_opt,~,exitflag] = fmincon(obective_func ,x0,[],[],[],[],lb,ub,...
        @(x) calc2_norm(x, row_indices, column_indices,dim,beta,issym,is_sparse,use_radius,use_my_svd),options);
    %         [],options);
    
   
    
    % ==== Create the output symmetric sparse matrix ====
    W = sparse(row_indices, column_indices,x_opt, dim,dim); 
    if ~is_sparse
        W = full(W);
    end
    
    if issym
        Wdiag = diag(diag(W));
        W = W + W' - Wdiag;
    end
    
    
        
    if use_radius && ~issym
        fprintf('start input: \t||A||_2 = %g\t ||A - 0||_F^2 = %g\n',  top2eigs(1)  , norm(A,'fro')  );
        fprintf('final output: \t R(W) = %g\t ||A - W||_F^2 = %g\n',  eigs(W,1)  , norm(W-A,'fro')  );
    else
        spect_norm = normest(W);
        fprintf('start input: \t||A||_2 = %g\t ||A - 0||_F^2 = %g\n',  start_2norm  , norm(A,'fro')  );
        fprintf('final output: \t||W||_2 = %g\t ||A - W||_F^2 = %g\n',  spect_norm  , norm(W-A,'fro')  );
    end

end

function [row_indices, column_indices, values] = extract_upper_half_inds(A)
  % This function takes the indices of the non zero element from the upper half
  % In case the matrix is symmetric you can we can use only the upper half 
  % of the matrix.
  %
    [row_indices, column_indices, values] = find(A);
    upper_half_ind = row_indices <= column_indices;
    row_indices = row_indices(upper_half_ind);
    column_indices = column_indices(upper_half_ind);
    values = values(upper_half_ind);
end

function [f,g] = func_objective(x,values)
    f = 0.5*norm(x - values).^2;
    g = x - values;
end

function [f,g] = func_objective_L1(x,values)
%     f = sum(values -x);
%     f = -sum(x); %sum(values) is const and x_ij < value_ij
%     g = -ones(size(values));
    
    f = sum(abs(x-values));
    g = sign(x - values);
end

function [c_leq,c_eq, g_leq, g_eq] = calc2_norm(X_vec, row_indices, col_indices,dim, beta, issym,is_sparse,use_radius,use_my_svd)
    % Computes the spectral constraint and the its gradient
    % We compute it only at the non-sparse locations
    %
    
%     if is_sparse
%         X_mat = sparse(non_zero_x, non_zero_y,X_vec, dim,dim);
%     else
%         X_mat = reshape(X_vec,dim,dim);
%     end
    
    X_mat = sparse(row_indices, col_indices,X_vec, dim,dim);
    if ~is_sparse
        X_mat = full(X_mat);
    end
    
    
    
    if issym
        sym_matrix = X_mat + X_mat' - diag(diag(X_mat));
%         sym_matrix = X_mat + triu(X_mat,1)';
        opts.issym = 1;
        opts.isreal = 1;
        [U,S] = eigs(sym_matrix,1,'lm', opts);
        V = U;
        S = abs(S);
%         [U,S,V] = svds(sym_matrix,1);

        % ========================
        % because I use only half the values of the matrix each value
        % contributes twice as much.
        % each value on the diagonal only contributes once
        symm_mult = 2*ones(size(row_indices));
%         find(non_zero_x==non_zero_y)
        symm_mult(row_indices==col_indices) = 1;
        g_leq = symm_mult.*U(row_indices,:).*V(col_indices,:);
%         g_leq = 2.*U(non_zero_x,:).*V(non_zero_y,:);

        
    else
        
        if use_radius
%             [U,S,V] = power_lowrank_vec( X_mat, rand(size(X_mat,2),1), rand(size(X_mat,1),1),...
%                     1000, 1e-10 ,issym,'k_eigs',true);
            
            [V,S] = eigs(X_mat,1,'lm');
            [U,S] = eigs(X_mat',1,'lm');
            
            sign_diag = sign(S);
            V = sign_diag*V;
            S = abs(S);
            g_leq = U(row_indices,:).*V(col_indices,:);
            g_leq = g_leq / dot(U,V);
        else
            if use_my_svd
                [U,S,V] = power_lowrank_vec( X_mat, rand(size(X_mat,2),1), rand(size(X_mat,1),1),...
                    1000, 1e-10 ,issym);
            else
                [U,S,V] = svds(X_mat,1);
            end
            g_leq = U(row_indices,:).*V(col_indices,:);
        end
        
    end
      
    c_leq = S - beta;
    c_eq = [ ];
    g_eq = [];
end
