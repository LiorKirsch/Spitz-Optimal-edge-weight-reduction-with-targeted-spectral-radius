function [A_projected,u,v] = project_on_bound_max_singular_val(A, max_val, is_sym,k) %,u_init,v_init)

%     max_iter_power_method = 3000;
% 
%     if ~exist('u_init','var')
%         u_init = rand(size(A,1),1);
%     else
%         k = size(u_init,2);
%     end
%     if ~exist('v_init','var')
%         v_init = rand(size(A,2),1);
%     end
%     [ u,s,v, it_num ] = power_method_svd( A, v_init, u_init, max_iter_power_method, 1e-16 ,is_sym);

% ======= uses svd for the first k

if ~exist('k','var')
    k = 1;
%     k = 'full';
end

if ismember(k, 'full');
    [ u,s,v ] = svd( A);
else
    [ u,s,v ] = svds( A,k);
end
         
    s_new = min( max_val , abs(s)) .*sign(s); 
    A_projected = A - u*(s - s_new)*v';

%     ==== Same as :  ====
%     A_projected = u*s_new*v';

    

    if is_sym
        A_projected = (A_projected+A_projected')/2; % Ensure symmetry.
    end
end 