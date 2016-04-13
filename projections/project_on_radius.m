function A_new = project_on_radius(A, beta)

[U,T] = schur(A,'complex');
% [U,T] = schur(A,'real');
eigen_vals = diag(T);
T_without_diag = T - diag(eigen_vals);

%==========================================
% change every eigne value that has a modulo that is bigger than beta 
% to have a module which is beta.
%==========================================

inds = ( beta < abs(eigen_vals));

above_eigen_vals = eigen_vals(inds);

% % eigen_vals( inds) = sign(eigen_vals(inds)) *beta;
% angles = angle(eigen_vals(inds));
% real_part = cos(angles) *beta;
% img_part  = sin(angles) *beta;
% eigen_vals( inds) = complex(real_part, img_part);

above_eigen_vals = above_eigen_vals ./ abs(above_eigen_vals) .*beta ;
eigen_vals( inds) = above_eigen_vals;


% T(inds,:) = bsxfun(@rdivide, T(inds,:) ,eigen_vals(inds)) *beta;


T = T_without_diag + diag(eigen_vals);

A_new = U*T*U';
A_new = real(A_new);

end