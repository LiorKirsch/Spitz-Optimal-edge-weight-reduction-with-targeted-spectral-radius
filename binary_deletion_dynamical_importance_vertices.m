function [H, df,remove_order] = binary_deletion_dynamical_importance_vertices(G,b,options)
% Greedy removal of vertices, one vertix at a time.
% input:
%   G - sparse undirected matrix
%   b - beta
%   options:
%      options.alg_type - 'first_order' or 'second_order' approximations
%      options.timelimit - the maximum number of seconds to run the method. 
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
    
    start_time = tic;

    is_sym = issymmetric(G);
    
    [d,i] = get_spectral_sorting(G,is_sym,alg_type);
    df = 0;
    remove_order = inf(size(G,1),2); % row,col,order
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
                string_to_print = print_output(string_to_print, abs(d), b, k);
            end
            if is_sym
                k = k+1;
                df = df + sum(G(i,:)); 
                df = df + sum(G(:,i)); 

                G(i,:) = 0; % remove incoming edges
                G(:,i) = 0; % remove outgoing edges

                remove_order(k,1)= i;
                remove_order(k,2)= k;

            else
                k = k+1;
                df = df + sum(G(i,:)); 
                df = df + sum(G(:,i)); 

                G(i,:) = 0; % remove incoming edges
                G(:,i) = 0; % remove outgoing edges

                remove_order(k,1)= i;
                remove_order(k,2)= k;
            end

            [d,i] = get_spectral_sorting(G,is_sym,alg_type); % recalculate the vertices order
        end
    end
    H = G;
    
    print_output(string_to_print, abs(d), b, k);
    fprintf('\n');
end

function string_to_print = print_output(former_string, radius, beta, k)
    fprintf(repmat('\b',1,length(former_string)));
    string_to_print = sprintf('%g ,(%d vertices)', radius/beta,k );
    fprintf('%s',string_to_print);
end

function [d,i,j] = get_spectral_sorting(G,is_sym,alg_type)
% We rank the vertices by computing the first or second order approximations
% of vertex removal
%
% ---- note ----
% The gradient should also be normalized by v*u in the
% case of non-symetric graphs but this has no effect on the ordering.
%
    
    switch alg_type
        case 'first_order'
            if  is_sym
                [u,d] = eigs(G,1); 
                ranking = u.*u;
            else
                [u,d] = eigs(G,1);
                [v,~] = eigs(G',1);
                ranking = abs(v.*u);
            end
        case 'second_order'
%             error('second_order not supported yet')
            if  is_sym
                [u,d] = eigs(G,1); 
                edge_sum = sum(G.*G',2);
                ranking = (-1 + edge_sum/ d^2) .*u.*u;
                ranking = abs(ranking); % is this nessecery ?
            else
                [u,d] = eigs(G,1);
                [v,~] = eigs(G',1);
                
                edge_sum = sum(G.*G',2);
                ranking = (-1 + edge_sum/ d^2) .*u.*v;
                ranking = abs(ranking);% is this nessecery ?
            end
        otherwise
            error('unkown parameter %s', alg_type);
    end
    
    [~, max_ind] = max(ranking);
    i = max_ind;
end

