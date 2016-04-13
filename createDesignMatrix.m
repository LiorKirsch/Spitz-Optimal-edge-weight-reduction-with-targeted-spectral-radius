function [D,edge_vals] = createDesignMatrix(A)
    num_edges = nnz(A);
    num_vertices = size(A,1);
    [rows,cols,edge_vals] = find(A);
    
    designs_rows = cat(1,1:num_edges,1:num_edges);
    designs_rows = designs_rows(:);
    designs_cols = cat(2,rows,cols)';
    designs_cols = designs_cols(:);
    
    D = sparse(designs_rows, designs_cols, ones(size(designs_cols)), num_edges, num_vertices);
    
end