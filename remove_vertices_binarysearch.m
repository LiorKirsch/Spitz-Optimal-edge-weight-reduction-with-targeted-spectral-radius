function modfied_A = remove_vertices_binarysearch( A, vertices_scores ,bound,tol,use_radius)
% Removing whole edges from A until ||A||_2 < bound.
% The edges are removed from big to small according to the size  A-W.
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

is_sym = isequal(A,A');
dim = size(A,1);

% [x_ind, y_ind, A_vals] = find(A);

% relevent_inds = sub2ind(size(A), x_ind, y_ind);
% W_vals = full(W(relevent_inds));

% [sorted_W_vals, sort_ind] = sort(W_vals);
[sorted_diff_vals, sort_ind] = sort(vertices_scores,'ascend');

ind_upper = 1;
ind_lower = dim;

norm2_upper = inf;
norm2_lower = 0;

fprintf('|u_ind - l_ind|= ');
string_to_print = '';
while (ind_lower -ind_upper ) >1
    
   ind_middle = round( (ind_lower + ind_upper)/2);
   
   modfied_A = remove_vertices_above_ind(A, ind_middle, sort_ind);
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


modfied_A = remove_vertices_above_ind(A, ind_lower, sort_ind);
% removed_edges = [x_ind(sort_ind(1:ind_upper)) ,y_ind(sort_ind(1:ind_upper))];

% if is_sym
%     % Insure that modified_A is symmetric by removing the symmetric edges
%     set_diff = setdiff(removed_edges,removed_edges(:,[2,1]),'rows');
%     edges_to_zero = set_diff(:,[2,1]);
%     set_diff = setdiff(removed_edges(:,[2,1]),removed_edges,'rows');
%     edges_to_zero = [edges_to_zero ; set_diff];
%     
%     modfied_A(edges_to_zero(:,1), edges_to_zero(:,2) ) = 0;
%     removed_edges = [removed_edges;edges_to_zero ];
% end


end

function modfied_A = remove_vertices_above_ind(A, ind, sort_ind)
%      indices_to_leave = sort_ind(ind:end);
     indices_to_remove = sort_ind(1:(ind-1));
     modfied_A = A;
     modfied_A(indices_to_remove,:) = 0;
     modfied_A(:,indices_to_remove) = 0;
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