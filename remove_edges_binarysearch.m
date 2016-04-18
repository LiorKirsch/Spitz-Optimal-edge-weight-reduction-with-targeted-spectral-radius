function [modfied_A,removed_edges] = remove_edges_binarysearch( A, W ,bound,tol,use_radius)
% Using a ranking to remove whole edges from A until ||A||_2 < bound.
% The edges are first ranked from big to small according to the size  A-W.
% Than, we remove the edges using binarysearch until find the smallest set
% which obeys the constraint ||A||_2 < bound.
%
% Input:
%      A     - the original matrix.
%      W     - a submatrix whose values are bounded by A and its 2norm by
%              bound.
%      bound - a bound for the 2norm
%      tol   - a tolerance for the calculation of 2norm at each step.
%              we use this tolerance and we increase it as we get closer to
%              the target point.
%      use_radius - should we use the spectral-radius or the spectral-norm

if ~exist('tol','var')
    tol = 1e-10;
end
if ~exist('use_radius','var')
    use_radius = true;
end

is_sym = issymmetric(A);
dim = size(A,1);
total_edges = nnz(A);

[x_ind, y_ind, A_vals] = find(A);

relevent_inds = sub2ind(size(A), x_ind, y_ind);
W_vals = full(W(relevent_inds));

% [sorted_W_vals, sort_ind] = sort(W_vals);
[sorted_diff_vals, sort_ind] = sort(A_vals - W_vals,'descend');

ind_upper = 1;
ind_lower = length(W_vals);

norm2_upper = inf;
norm2_lower = 0;

fprintf('|u_ind - l_ind|= ');
string_to_print = '';
while (ind_lower -ind_upper ) >1
    
   ind_middle = round( (ind_lower + ind_upper)/2);
   
   modfied_A = remove_edges_above_ind(x_ind, y_ind, A_vals, dim, ind_middle, sort_ind);
   [is_below, norm2_modA] = check_for_A_below_bound(modfied_A,bound,tol,is_sym,use_radius);
   
   if is_below
       ind_lower = ind_middle;
       norm2_lower = norm2_modA;
   else
       ind_upper = ind_middle;
       norm2_upper = norm2_modA;
   end
   
   if abs(norm2_upper - norm2_lower)*1e-16  < tol 
       tol = tol *0.1;
   end

   fprintf(repmat('\b',1,length(string_to_print)));
   string_to_print = sprintf('%d ', ind_lower -ind_upper );
   fprintf('%s',string_to_print);
end


modfied_A = remove_edges_above_ind(x_ind, y_ind, A_vals, dim, ind_lower, sort_ind);
removed_edges = [x_ind(sort_ind(1:ind_upper)) ,y_ind(sort_ind(1:ind_upper))];

if is_sym
    % Insure that modified_A is symmetric by removing the symmetric edges
    edges_to_zero = setdiff(removed_edges(:,[2,1]),removed_edges,'rows');
    
    modfied_A(edges_to_zero(:,1), edges_to_zero(:,2) ) = 0;
    removed_edges = [removed_edges;edges_to_zero ];
end

num_edges_removed = size(removed_edges,1);
fprintf('\t ,(%d edges, %g%%)\n', num_edges_removed, num_edges_removed/total_edges);

end

function modfied_A = remove_edges_above_ind(x_ind, y_ind, A_vals, dim, ind, sort_ind)
     indices_to_leave = sort_ind(ind:end);
     modfied_A = sparse(x_ind(indices_to_leave), y_ind(indices_to_leave), A_vals(indices_to_leave), dim,dim);
end
function [is_below, A_2norm] = check_for_A_below_bound(A,bound,tol,is_sym,use_radius)
%     max_power_iter = 1000;
%     [~,A_2norm,~,~] = power_lowrank_vec( A, rand(size(A,2),1), rand(size(A,1),1),...
%          max_power_iter, tol ,is_sym, 'norm',use_radius);
    if use_radius
        A_2norm = abs(eigs(A,1));
    else
        A_2norm = normest(A);
    end
    is_below = A_2norm < bound;   
end