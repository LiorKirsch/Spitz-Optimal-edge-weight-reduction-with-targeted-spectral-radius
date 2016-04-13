function [W,vertices_weight,exitflag] = near_bounded_interior_p_vertices(A, beta, P, options)
%     This function computes the nearest bounded matrix to the matrix A.
%     Each element in the output matrix W is bounded above by the value of A
%     and bounded bellow by 0.
%     The largest eigenvalue of W is bounded by beta.
%     If the matrix A is symmetric we only use half the variables. 
%
%     Input: 
%           A  - sparse adjacancy matrix.
%           beta - the bound for the first eigen value (the spectral norm)
%           P - a matrix which controls the effect of pairs of vertex pair on each edge
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
%
%     Output: 
%           W - the nearest matrix
%           vertices_weights - the weights assigned to each of the vertices
%           exitflag - the end status of the optimization.
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
%           W = near_bounded_interior_p_vertices(A, beta, []);
%
if exist('P','var')
    if isempty(P)
        P = createVerticesInfluenceMatrix(A,'equal_weight');
    end
else
    P = createVerticesInfluenceMatrix(A,'equal_weight');
end
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

% find vertices which does affect the objective
    
    
    % ==== Perform some checks and statistics ====
    [dim, dim2] = size(A);
    assert(dim==dim2,'A should be a squared matrix');
    
    is_sparse = issparse(A);
    issym = issymmetric(A); % not sure this is supported in the new format with the P matrix.
    assert(issym, 'not sure this is supported in the new format with the P matrix.');
    issym = false ; 
    
    if issym
        fprintf('using symmetric properties\n');
    end
    opts.issym = issym;

    top2eigs = eigs(A,2,'lm', opts);
    start_2norm = normest(A);
    fprintf('Spectral radius of the original problem is %g\n', abs(top2eigs(1)));
    fprintf('Spectral  norm  of the original problem is %g\n', start_2norm );
    fprintf('Spectral gap is  %g\n', ( abs(top2eigs(1)) - abs(top2eigs(2)) ));

    if issym
        [row_indices, column_indices, values] = extract_upper_half_inds(A);
    else
        [row_indices, column_indices, values] = find(A);
    end
   
    % ===   using the fact that ||A||_2 <= ||A||_F or that R(A) <= ||A||_F ===
    x0 = values * min(1, beta/(norm(values) + eps));
    values_P = P(sub2ind(size(P), row_indices, column_indices));
    D = get_matrix_values(dim,row_indices, column_indices,values, values_P);
    x0 = D'*x0;
    x0 = zeros(size(x0));

    
    % ====use only vertices which have an affect====
    vertices_which_affect = full(sum(D,1) > 0);
    x0 = x0(vertices_which_affect);
    D = D(:,vertices_which_affect);
    
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

    lb = zeros(dim,1);
    ub = ones(dim,1);
    ineq_A = [];
    ineq_b = [];      
 
    
    switch objective
        case 'L1'
            obective_func = @(x) func_objective_L1(x,values,D);
        case 'L2'
            obective_func = @(x) func_objective(x,values,D);
        otherwise
            error('unkown objective - %s, should be ''L1'' or ''L1'' ', objective);
    end
   
    
     %%%%% yonatan added new line here because x0 was not fisible point
     %%%%% since the inequility Wi+Wj < Aij was not satisfied
%     tr = 0.5*max(max(A)); 
%     x0(x0>tr)=tr;
    %%%%%%
    %%%%%% 
    
    [x_opt,~,exitflag] = fmincon(obective_func ,x0,ineq_A,ineq_b,[],[],lb,ub,...
        @(x) calc2_norm(x, D,row_indices, column_indices,dim,beta,issym,is_sparse,use_radius,use_my_svd),options);
    %         [],options);
    
   
    
    % ==== Create the output matrix ====
    vertices_weight = nan(length(vertices_which_affect),1);
    vertices_weight(vertices_which_affect) = x_opt;
    
    x_opt = D*x_opt;
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

function [f,g] = func_objective(x,values,D)
    x_edge = D*x;
    
    f = 0.5*norm(x_edge - values).^2;
    g = D'*(D*x-values);
end

function [f,g] = func_objective_L1(x,values,D)
    x_edge = D*x;

    f = sum(abs(x_edge - values));
    g = D'*sign(x_edge - values);
end

function [c_leq,c_eq, g_leq, g_eq] = calc2_norm(X_vec, D,row_indices, col_indices,dim, beta, issym,is_sparse,use_radius,use_my_svd)
    % Computes the spectral constraint and the its gradient
    % We compute it only at the non-sparse locations
    %
    

    X_vec = D*X_vec;
    
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
            
%             
            [U,S] = eigs(X_mat,1,'lm');
            V = U;
%             [V,S] = eigs(X_mat',1,'lm');

            sign_diag = sign(S);
            V = sign_diag*V;
            S = abs(S);
            
        else
            if use_my_svd
                [U,S,V] = power_lowrank_vec( X_mat, rand(size(X_mat,2),1), rand(size(X_mat,1),1),...
                    1000, 1e-10 ,issym);
            else
                [U,S,V] = svds(X_mat,1);
            end
        end
        g_leq = U(row_indices,:).*V(col_indices,:);
    end
      
    c_leq = S - beta;
    c_eq = [ ];
    g_eq = [];
    
    g_leq = D'*g_leq;  % the inner function gradient
end

function D = get_matrix_values(num_vertices,row_indices, column_indices,values, p_vals)
    num_edges = length(row_indices);
    designs_rows = cat(1,1:num_edges,1:num_edges);
    designs_rows = designs_rows(:);
    designs_cols = cat(2,row_indices,column_indices)';
    designs_cols = designs_cols(:);
    
    values_cols = cat(2,values,values)';
    values_cols = values_cols(:);
    
    values_p = cat(2,p_vals,1 - p_vals)';
    values_p = values_p(:);
    
    values_cols = values_cols .* values_p;
    
    D = sparse(designs_rows, designs_cols, values_cols, num_edges, num_vertices);
end
