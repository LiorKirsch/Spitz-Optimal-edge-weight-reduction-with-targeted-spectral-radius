function [H, df,remove_order] = binary_deletion_dynamical_importance(G,b,options)
% Greedy removal of edges, one edge at a time.
% input:
%   G - sparse undirected matrix
%   b - beta
%   options:
%      options.alg_type - 'first_order' or 'second_order' approximations
%      options.use_cached_eigvecs - binary, whether to use cached vectors for the
%            eigen decomposition.
%      options.timelimit - - the maximum number of seconds to run the method. 
%              When exceeded return NaN.
%
% output:
%   H - the graph after the edges are removed
%   df - the total sum of the removed edges
%

    if ~exist('options','var')
        options = struct();
    end
    if isfield(options,'time_limit')
        time_limit = options.time_limit;
    else
        time_limit = inf;
    end
    if isfield(options,'alg_type')
        alg_type = options.alg_type;
    else
        alg_type = 'first_order';
    end
    if isfield(options,'use_cached_eigen_vec')
        use_cached_eigen_vec = options.use_cached_eigen_vec;
    else
        use_cached_eigen_vec = false;
    end
    

    start_time = tic;
    is_sym = issymmetric(G);
    [d,i,j,u,v] = get_spectral_sorting(G,is_sym,alg_type,false); % first time don't use cached vecs
    df = 0;
    total_edges = nnz(G);
    remove_order = inf(total_edges,3); % row,col,order
    k=0;
    fprintf('%s dynamical importance R(G)/beta: ',alg_type);
    string_to_print = '';
    while abs(d) > b
        if time_limit < toc(start_time) 
            disp('======= time limit reached, stopping! =======');
            G = nan;
            df = nan;
            remove_order = nan;
            break; 
        else
            if ~mod(k,10) % display progress towards beta
                string_to_print = print_output(string_to_print, abs(d), b, k, total_edges);
            end
            if is_sym
                k = k+2;
                df = df + 2*G(i,j); %*2 because g is symmetric

                G(i,j) = 0; % delete the edge
                G(j,i) = 0; % delete also the symmetric edge

                remove_order(k-1,1)= i;
                remove_order(k-1,2)= j;
                remove_order(k-1,3)= k/2;
                remove_order(k,1)= j;
                remove_order(k,2)= i;
                remove_order(k,3)= k/2;

            else
                k = k+1;
                df = df + G(i,j); %*2 because g is symmetric

                G(i,j) = 0; % delete the edge

                remove_order(k,1)= i;
                remove_order(k,2)= j;
                remove_order(k,3)= k;
            end

            % recalculate the edge order
            [d,i,j,u,v] = get_spectral_sorting(G,is_sym,alg_type,use_cached_eigen_vec,u,v); 
        end
    end
    H = G;
    
    print_output(string_to_print, abs(d), b, k, total_edges);
    fprintf('\n');
end

function string_to_print = print_output(former_string, radius, beta, k, total_edges)
    fprintf(repmat('\b',1,length(former_string)));
    string_to_print = sprintf('%g ,(%d edges , %g%%)', radius/beta,k ,k/total_edges);
    fprintf('%s',string_to_print);
end

function [d,i,j,u,v] = get_spectral_sorting(G,is_sym,alg_type,use_cached_eigvecs,old_u,old_v)
% We rank the edges by computing the first or second order approximations
% of edge removal
%
% ---- note ----
% The gradient should also be normalized by v*u in the
% case of non-symetric graphs but this has no effect on the ordering.
%
    [row,col,vals] = find(G);
    power_method_tol = 1e-16;
    power_method_max_iter = 500;
    
    switch alg_type
        case 'first_order'
            if  is_sym
                if use_cached_eigvecs
                    [ u, d ] = power_method_k(G, 1, old_v, power_method_max_iter, power_method_tol);
                else
                    [u,d] = eigs(G,1); 
                end
                v = u;
                
                ranking = abs(u(row).*u(col).*vals);
            else
                if use_cached_eigvecs
                    [ u, d ] = power_method_k(G, 1, old_u, power_method_max_iter, power_method_tol);
                    [ v, ~ ] = power_method_k(G', 1, old_v, power_method_max_iter, power_method_tol);
                else
                    [u,d] = eigs(G,1);
                    [v,~] = eigs(G',1);
                end
                ranking = abs(v(row).*u(col).*vals);
            end
        case 'second_order'
            if  is_sym
                [u,d] = eigs(G,1);
                v = u;
                ranking = 2*vals .*u(row).*u(col)/d + ...
                          vals.^2 .* ( u(row).^2 + u(col).^2) /d^2;
            else
                [u,d] = eigs(G,1);
                [v,~] = eigs(G',1);
                only_self_loops = row==col;

                ranking = vals.*v(row).*u(col)/d ;

                % === less readable more efficent ===
%                 row_loops = row(only_self_loops);
%                 col_loops = col(only_self_loops);
%                 second_order_fix = zeros(size(ranking));
%                 second_order_fix(only_self_loops) = vals(only_self_loops).^2 .*v(row_loops).*u(col_loops)/d^2;
                
                % === more readable less efficent ===
                second_order_fix = vals.^2 .*v(row).*u(col)/d^2;
                second_order_fix(~only_self_loops) = 0;
                
                ranking = ranking + second_order_fix;
            end
        case 'k_approximation'
            G = full(G);
            if use_cached_eigvecs
                [ u, d ] = power_method_k(G, 1, old_v, power_method_max_iter, power_method_tol);
            else
                [u,d] = eigs(G,1); 
            end
            v = u;
            %compute walk length k
            n = size(G,1);
            k=2*log(n)-2;
            iter=ceil(log2(k));

            if k<2
                k=2;
                iter=1;
            end
            
            %Compute matrix adm^k
            adm = G;
            admk=adm;
            for mult=1:iter
                admk=admk*admk;
                min_pos_val=min(nonzeros(admk));
%                 admk=admk/min_pos_val;
            end
            admk=admk.*adm; % multiple the number of k-walks with the cost of each edge

            ranking = admk(sub2ind(size(adm),row,col));

        otherwise
            error('unkown parameter %s', alg_type);
    end
    
    [~, max_ind] = max(ranking);
    i = row(max_ind);
    j = col(max_ind);
end

