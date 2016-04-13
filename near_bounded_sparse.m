function [Y,iter] = near_bounded_sparse(A,eigen_bound,options)
%   [Y,iter] = near_bounded_sparse(A,eigen_bound,options)
%   finds the nearest reduction with bounded specturn matrix.
%   Input: 
%       A - a sparse matrix. 
%       eigen_bound  - an upper bound for the largest eigen value.
%       options.TolFun  - is a convergence tolerance, which defaults to 16*EPS. 
%       options.MaxIter - the maximum number of iterations (default 100, but may
%           need to be increased, depends on the number of eigen values that
%           we need to cap).
%       options.Display - a flag which contronls the print (0 - no print, 1 - print,
%           2 - print also includes the spectral norm.
%       options.InteriorDist - since dykstra can go on and on and only converge
%           after infinit number of steps, we use a small dist tp find a point 
%           inside the convex group of matrices with bounded spectrum
%           (default 1e-10)
%       options.UseRadius - if A is not symmetric we can use the svd (default) or
%           the schur factorization to get the bounds on the radius
%           true - specifices using schur, false (default) uses svd
%       options.K    - the number of eigen vectors to use for the projection
%           (default 1)
%
%   Output: 
%       Y - the nearest matrix
%       iter - the number of iterations used.
%
%  We use alternating projections with dykstra correction to find the closest point.
%  We first remove the zero columns and zero rows. Then, we use alternating
%  projection between the bounds projection and the spectral bound
%  projection. 
%
%   Author:
%      Lior Kirsch 2015

% ------- Check parameters
if ~exist('options','var') options = struct; end
if ~isfield(options,'Display') options.Display='final';end
if ~isfield(options,'MaxIter') options.MaxIter=3000;end
if ~isfield(options,'TolFun') options.TolFun=1e-6;end
if ~isfield(options,'UseRadius') options.UseRadius=false;end
if ~isfield(options,'InteriorDist') options.InteriorDist=1e-10;end
if ~isfield(options,'K') options.K=1;end
if ~isfield(options,'use_cached_eigen_vec') options.use_cached_eigen_vec=true;end
if ~isfield(options,'power_method_tol') options.power_method_tol=1e-16;end

assert( issparse(A) ,'A must be sparse');
warning('off','MATLAB:normest:notconverge');
string_to_print = '';

max_power_iter = 500;
is_sym = issymmetric(A);
if is_sym && (ismember(lower(options.Display)  ,{'iter','final','details'}))
    fprintf('using symmetric properties\n');
end

%======== remove zero column and zero rows =============
[A, nonzero_rows, nonzero_col] = get_zero_row_columns(A);

[row_indices, column_indices, bound_values] = find(A);
assert( all(bound_values >=0) ,'all entries in A must greater or equal to 0');
Y_values = bound_values;
dS2 = 0; dS_values = 0; dS_values_old = 0;

X = A; Y = A;
iter = 0;
rel_diffX = inf; rel_diffY = inf; rel_diffXY = inf; rel_diffXY2 = inf;

% options.K =3;
% singular_vals = svds(A,options.K);
% options.K = max(1, sum(singular_vals > eigen_bound) );

s = zeros(options.K); s_new = zeros(options.K);
if is_sym
    [u,~] = eigs(A, options.K);
    v = u;
else
    [u,~] = eigs(A, options.K);
    [v,~] = eigs(A', options.K);
end
% ====== since s is zeros this inits the slack variable L with 0.
% ====== the init of u,v is just for the case when cache is used.
u_old = u;
v_old = v;
s_old = s;
s_new_old = s_new;

if ismember(lower(options.Display)  ,{'iter','final','details'})
  
  head_string = 'i:       ||Y - Yold||_2   ||X - Xold||_2   ||X - Y||';
  fprintf('%s',head_string);
  if ismember(lower(options.Display)  ,{'details'}) 
      fprintf('  ||Y||_2 - b');
      [~,Y_norm,~,~] = power_lowrank_vec( Y, u, v,...
         max_power_iter, min(1e-10,options.InteriorDist/100) ,is_sym, 'norm',options.UseRadius);
      fprintf('\n-:%s%9.2e', repmat(' ',1,length(head_string)),...
                   Y_norm - eigen_bound);   
  end
  fprintf(' (tolerance %g)\n',options.TolFun);
end


while max([rel_diffX rel_diffY rel_diffXY]) > options.TolFun && (iter <= options.MaxIter)
 
   u_old_old = u_old;
   v_old_old = v_old;
   s_old_old = s_old;
   s_new_old_old = s_new_old;
   u_old = u;
   v_old = v;
   s_old = s;
   s_new_old = s_new;
   dS_values_old_old = dS_values_old;
   dS_values_old = dS_values;
   Yold = Y;
   Y_values_old = Y_values;
   
   if options.use_cached_eigen_vec
       power_u0 = u; 
       power_v0 = v; 
   else
       power_u0 = rand(size(A,1),options.K); 
       power_v0 = rand(size(A,2),options.K); 
   end
   
   % This is the same as: [u,s,v] = svds(Y + u*(s-s_new)*v',options.K);
   % ... = proj1(Y + dS);
   [ u,s,v, it_num ] = power_lowrank_vec_sparse_low(Y ,u,s - s_new,v, power_u0,power_v0, max_power_iter, options.power_method_tol ,is_sym,options.UseRadius);
   
    
%    s_new = min( eigen_bound - options.InteriorDist , abs(s)) .*sign(s);
   
   s_diag = diag(s);
   s_new_diag = min( eigen_bound - options.InteriorDist , abs(s_diag));    
   s_new = s - diag(s_diag) + diag(s_new_diag);%in the case where s is triag and not diag    
   
%    dS is the dykstra correction.   
%    dS = (Y + dS) - X;
%    dS = u*(s - s_new)*v';

 
   % Here we compute the values of X just in the non-sparse locations
   % The values of the correction var dS2 are also computed just at the non
   % sparse locations.
   
   dS_values = get_x_vals(row_indices, column_indices,u,s -s_new,v);

   
   % ===== compute || X - X_old||_F ====
   rel_diffX = norm(Y_values - Y_values_old - dS_values + dS_values_old_old)^2;
   rel_diffX = rel_diffX + norm((s - s_new) - (u'*u_old_old) *(s_old_old - s_new_old_old)*(v_old_old'*v),'fro')^2;
   rel_diffX = rel_diffX - norm(dS_values - dS_values_old_old)^2;
   rel_diffX = sqrt(rel_diffX);
   rel_diffX = rel_diffX / sqrt(norm(Y_values - dS_values)^2 + norm(s - s_new)^2 -norm(dS_values)^2 );
   % ===================================

   
   % X_vals =    (Input to proj_svd)            - values_out_of_projections
   X_values = (Y_values + dS_values_old) - dS_values;
   Y_values = project_on_range(X_values + dS2, bound_values);
      
   dS2_old = dS2;
   dS2 = X_values - Y_values  + dS2;
   Y = sparse(row_indices, column_indices,Y_values, size(A,1), size(A,2));
   
   rel_diffY = norm(Y-Yold,'fro')/norm(Y,'fro');
   rel_diffXY = norm((s - s_new) - (u'*u_old) *(s_old - s_new_old)*(v_old'*v),'fro')/norm((s - s_new),'fro');
   
%    robust_crit = norm(s-s_new)^2 - norm(s_old-s_new_old)^2 +  norm(dS2)^2 - norm(dS2_old)^2;
   if ismember(lower(options.Display)  ,{'iter','details'})
      fprintf(repmat('\b',1,length(string_to_print)));
      string_to_print = sprintf('%2.0f:\t  %3.2e\t    %3.2e\t   %3.2e \t %3.2e', ...
           iter, rel_diffY, rel_diffX, rel_diffXY);%, robust_crit);
      if ismember(lower(options.Display)  ,{'details'})
          [~,Y_norm,~,~] = power_method_svd( Y, rand(size(Y,2),1), rand(size(Y,1),1),...
            max_power_iter, min(1e-10,options.InteriorDist/100) ,is_sym, 'norm',options.UseRadius);
          string_to_print = sprintf('%s  %3.2e', string_to_print,Y_norm - eigen_bound);
      end
      fprintf('%s',string_to_print);
   end
   iter = iter + 1;
end

if iter > options.MaxIter
    fprintf('\nStopped after %d iterations. Try increasing MAXITS.\n',options.MaxIter );
else
    fprintf('\nStopped after %d iterations with bound dist %g. \n',iter,options.InteriorDist );
end

 
%======== add back the  zero columns and rows =============
    Y = add_zeros(Y,nonzero_rows,nonzero_col);
    warning('on','MATLAB:normest:notconverge');
end


function [clean_A, nonzero_rows, nonzero_col] = get_zero_row_columns(A)
    nonzero_rows = any(A,2);
    nonzero_col = any(A,1);
    
    nonzero_col_or_row = nonzero_rows | nonzero_col';
    nonzero_rows = nonzero_col_or_row;
    nonzero_col = nonzero_col_or_row;
    clean_A = A(nonzero_col_or_row, nonzero_col_or_row);
end
function new_A = add_zeros(A,nonzero_filter_rows,nonzero_filter_cols)
% Adds back the zero columns and rows that were removed at the begining
%

%     new_A = sparse(length(nonzero_filter_rows),length(nonzero_filter_cols) );
%     new_A2(nonzero_filter_rows, nonzero_filter_cols) = A;
    
    translation_row = find(nonzero_filter_rows);
    translation_col = find(nonzero_filter_cols);
    [A_row, A_col, A_val] = find(A);
    new_A_row = translation_row(A_row);
    new_A_col = translation_col(A_col);
    new_A = sparse(new_A_row, new_A_col, A_val,...
        length(nonzero_filter_rows),length(nonzero_filter_cols));
end

function X_values = get_x_vals(row_indices, column_indices,u,s,v)
% Returns the values in the indices (row_indcies(i),column_indcies(i))
% for the matrix u*s*v'.
% It does this without explicitly computing the matrix which might be too
% large.
%

    % out_values = arrayfun(@(x,y) u(x)*u(y), row_indices, column_indices);

    % ============== This is actually:  x = u * s *v; and then take ==========
    % ============== only the values using the filter               ==========
    X_values = u(row_indices,:).*v(column_indices,:);
    X_values = X_values * diag(s);
    
    % if length(s) ==1
    %     X_values = u(row_indices,1).*v(column_indices,1) * s;
    % else
    %     addition_matrix = u* s *v';
    %     inds = sub2ind(size(addition_matrix),row_indices,column_indices);
    %     X_values = addition_matrix(inds);
    % end
end