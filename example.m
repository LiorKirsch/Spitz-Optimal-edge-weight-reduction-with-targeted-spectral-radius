  addpath('projections/'); 
  addpath('lbfgs/');
  
  % create a random undirected network
  dim = 1000;       num_nnz = dim^2 /100; 
  
  denisty = num_nnz/(dim^2);
  A = sprand(dim,dim, denisty);
  A = A - diag(diag(A));
  A = (A + A')/2;

% ==== Some optimization parameters ====
options = optimset();
options.Display = 'iter' ;%'final'; %'iter';
options.MaxIter = 3000;
options.TolFun = 1e-6;
options.InteriorDist = 1e-5;
options.K = 1;
options.Objective = 'L2';

eval_func = @(x) normest(x);

beta  = round(eval_func(A) * 0.7);
fprintf('%d singular values above bound \n', sum(svds(A,100) > beta) );
fprintf('%d non-zero values in A\n', nnz(A));    

% solve the edges problem
X_dykstra = near_bounded_sparse(A,beta, options);
X_intp = near_bounded_interior_p(A, beta,options);
X_dynamical_importance = binary_deletion_dynamical_importance(A, beta,options);

% solve the vertices problem
P = createVerticesInfluenceMatrix(A,'equal_weight');
X_dykstra_v = near_bounded_sparse_vertices(A,beta, P, options);
X_intp_v = near_bounded_interior_p_vertices(A, beta, P, options);
X_dynamical_importance_v = binary_deletion_dynamical_importance_vertices(A, beta,options);


fprintf('\n====================== edges ==================================\n');
fprintf('start input       : \t||A||_2 = %g\t ||A||_F = %g\t b= %g\n',  eval_func(A) , norm(A,'fro') , beta);
fprintf('Spitz dykstra: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_dykstra)  , norm(X_dykstra-A,'fro') );
fprintf('Spitz interior-points: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_intp)  , norm(X_intp-A,'fro') );
fprintf('Dynamical Importance: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_dynamical_importance)  , norm(X_dynamical_importance-A,'fro') );

fprintf('\n====================== vertices ==================================\n');
fprintf('start input       : \t||A||_2 = %g\t ||A||_F = %g\t b= %g\n',  eval_func(A) , norm(A,'fro') , beta);
fprintf('Spitz dykstra: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_dykstra_v)  , norm(X_dykstra_v-A,'fro') );
fprintf('Spitz interior-points: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_intp_v)  , norm(X_intp_v-A,'fro') );
fprintf('Dynamical Importance: \t||W||_2 = %g\t ||A - W||_F = %g\n',  eval_func(X_dynamical_importance_v)  , norm(X_dynamical_importance_v-A,'fro') );



