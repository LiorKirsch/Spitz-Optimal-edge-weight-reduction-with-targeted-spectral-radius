function [ u,s,v, it_num ] = power_lowrank_vec( A, v_init, u_init, it_max, tol ,is_sym, stop_crit,use_radius)
% Find the largest singular values and singular vectors using the power method.
% Find the largest (in absolute) eigen values and eigen vectors using the power method.
%
% Input:
%      A - the matrix 
%      v_init - init vector for the right singular vector
%      u_init - init vector for the left singular vector
%      it_num - maximum number of iteration for the power method
%      tol - tolerance for the power method
%      is_sym - (default true) is the matrix A symetric
%      stop_crit - can be "norm" or "k_eigs" (default: "k_eigs")
%      use_radius - if A is not symmetric we can use the svd (default) or
%                   the schur factorization to get the bounds on the radius
%                   true - specifices using schur, false (default) uses svd

    if ~exist('is_sym','var')
        is_sym = true;
    end
    if ~exist('stop_crit', 'var')
      stop_crit = 'k_eigs';
    end
    if ~exist('use_radius', 'var')
      use_radius = false;
    end

    if is_sym
        k = size(u_init,2);
        [ v, s, it_num ] = power_method_k(A,k, v_init, it_max, tol,stop_crit );
        u = v;
    else
        k = size(u_init,2);
        if use_radius
            [ u, s, it_num ] = power_schur_k(A,k, v_init, it_max, tol,stop_crit );
            v = u;
        else
            [ u, s,v, it_num ] = power_svd_k(A,k, u_init,v_init, it_max, tol,stop_crit );
        end
    end
    
end    

