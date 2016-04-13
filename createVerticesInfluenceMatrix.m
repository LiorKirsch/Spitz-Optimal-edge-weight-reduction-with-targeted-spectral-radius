function P = createVerticesInfluenceMatrix(A,scenario)
% builds the vertices influence matrix for symetric matrices only
% Input:
%       A - an adjcancy matrix
%       scenario - how to set the contribution of vertices to the edges
%            'equal_weight','hubs_focused','hubs_focused_binary','small_degree_focused','small_degree_focused_binary'
% Output:
%       P - [n,n] matrix where P_ij is the influence that vertex i has on
%           edge Aij
            
    if ~exist('scenario','var')
        scenario = 'equal_weight';
    end
    fprintf('creating vertices-influence matrix using %s scenario\n', scenario);
    
    num_vertices = size(A,1);
    [rows,cols,values] = find(A);
    
    switch scenario
        case 'equal_weight'
            entries_vals = 0.5 * ones(size(rows));
        case 'hubs_focused'
            vertices_degree = sum(A,1);
            entries_vals = vertices_degree(rows) ./ (vertices_degree(rows) + vertices_degree(cols));
        case 'hubs_focused_binary'
            vertices_degree = sum(A,1);
            entries_vals = vertices_degree(rows) ./ (vertices_degree(rows) + vertices_degree(cols));
            entries_vals = round(entries_vals);
        case 'small_degree_focused'
            vertices_degree = sum(A,1);
            entries_vals = vertices_degree(cols) ./ (vertices_degree(rows) + vertices_degree(cols));
        case 'small_degree_focused_binary'
            vertices_degree = sum(A,1);
            entries_vals = vertices_degree(cols) ./ (vertices_degree(rows) + vertices_degree(cols));
            entries_vals = round(entries_vals);
        otherwise 
            error('unknown matrix creation scenario - %s',  scenario);
    end
 
    P = sparse(rows, cols, entries_vals, num_vertices, num_vertices);
end