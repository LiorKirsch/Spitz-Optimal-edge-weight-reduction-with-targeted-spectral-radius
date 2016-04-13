function [ u,s,v, it_num ] = power_lowrank_vec_sparse_low(A,u,s,v, power_u0,power_v0, it_max, tol ,is_sym,use_radius)
% Find the largest singular values and singular vectors using the power method.
% This function is speciatlly tuned for the case that A is the sum of a
% sparse matrix and a low rank matrix
%
% Input: Z = A + u*s*v (A is sparse and u,v are basis vector)
%      A - the sparse matrix  (m x n)
%      u - the left low rank basis matrix  (m x k)
%      s - the singular value matrix       (k x k)
%      v - the right low rank basis matrix (n x k)
%      power_v0 - init vector for the right singular vector
%      power_u0 - init vector for the left singular vector
%      it_max - maximum number of iteration for the power method
%      tol - tolerance for the power method
%      is_sym - (default true) is the matrix Z symetric
%      use_radius - if A is not symmetric we can use the svd (default) or
%                   the schur factorization to get the bounds on the radius
%                   true - specifices using schur, false (default) uses svd
%      
% We use a variant of the power method to compute rank k svd for Z.
% We exploit the structure of Z which is "sparse" + "low rank"
% In case the matrix Z is symetric the power method computation is 
%  much simpler since the left and right eigen values are the same.
%

it_num = nan;
   if ~exist('is_sym','var')
        is_sym = true;
   end
   if ~exist('use_radius','var')
        use_radius = false;
   end
       
   if is_sym

        k = size(u,2);
        
        [ v, s, it_num ] = power_method_k_sparse_low(A, u, s,v, k, power_v0, it_max, tol,is_sym,size(A,1) );
        u = v;
        
%       this is the same as:  
%       [u,s,v] = svds(A + u*s* v' ,k);

   else
        k = size(u,2);
         
        if use_radius
            [ v, s, it_num ] = power_schur_k_sparse_low(A, u, s,v, k, power_v0, it_max, tol,is_sym,size(A,1) );
            u = v;
        else
            [u, s,v, it_num ] = power_svd_k_sparse_low(A ,k, u, s,v,  power_u0,power_v0, it_max, tol );
%             [u,s,v] = svds(A + u*s* v' ,k);
        end
      
   end
   
   
   % ==== this just make the next step easier                   ====
   % ==== even if the eigenvalues are negative, flip their sign ====
   sign_diag = sign(diag(s))';
   v = bsxfun(@times,v, sign_diag);
   s = bsxfun(@times,s, sign_diag);

   
end