function [ u, s,v, it_num ] = power_svd_k_sparse_low(A, k, u_low, s_low, v_low, u0,v0, it_max, tol )

% Find the k largest singular values and singular vectors using the power method.
%   uses the fact that Z = "sparse" + "low rank"
%   we never compute Z explicitly always leave it as the sum of two
%   matrices.
%
% Input: 
%      A - the matrix  (m x n)
%      k - the number of eigen values to find
%      u_low - the left part of the low rank matrix    (m x k)
%      s_low - the singular value part of the low rank (k x k)
%      v_low - the right part of the low rank matrix   (n x k)
%      u0 - init vector for the left singular vector   (m x k)
%      v0 - init vector for the right singular vector  (n x k)
%      it_max - maximum number of iteration for the power method
%      tol - tolerance for the power method
%

  debug = 0;
  min_iter = 10;

  if ( debug )
    fprintf ( 1, '\n' );
    fprintf ( 1, '     IT      Lambda          Delta-Lambda    Delta-Y\n' );
    fprintf ( 1, '\n' );
  end
  
  assert( size(u0,2) == k, 'u0 should be of size n X k') ;
  assert( size(v0,2) == k, 'v0 should be of size n X k') ;
  s = 0;
  u = u0;
  v = v0;
  it_num = 0;
  val_dif = 0.0;
  A_t = A';
  v_low_t = v_low';
  u_low_t = u_low';
  
%   norm_A = norm(A + u_low*s_low*v_low','fro');
%   norm_A1 = norm(A + u_low*s_low*v_low',1);
  norm_A2 = calc_fro_norm(A, u_low, s_low, v_low);
  tol_adjusted = tol*norm_A2;

  if ( debug )
    fprintf ( 1, '  %5d  %14e  %14e\n', it_num, s, val_dif );
  end

  for it_num = 1 : it_max

    s_old = s;
     if k==1
      av = A * v + u_low* s_low *(v_low_t * v);
      s = norm(av);
      u = av / s;

      au = A_t * u + v_low * s_low *(u_low_t *u);
      s = norm(au);
      v = au / s;
    else
      av = A * v + u_low* s_low *(v_low_t * v);
      [Q,R] = qr(av,0);
      u = Q;
      
      au = A_t * u + v_low * s_low *(u_low_t *u);
      [Q,R] = qr(au,0);
      v = Q;
      s = diag(diag(R));
     end
  
   
    val_dif = norm(au - v*s,1);
    

    if (debug)
      fprintf('  %5d  %14e  %14e \n', it_num, s, val_dif);
    end

    if ( val_dif <= tol_adjusted && min_iter <= it_num)
      break
    end 

  end

   v = v * sign(s);
   s = abs(s);
end

function fro = calc_fro_norm(A, u, s, v)
    [row_indices, column_indices, vals] = find(A);
    low_vals_in_nonzero = (u(row_indices,:).*v(column_indices,:)) * diag(s);

    fro = norm( vals + low_vals_in_nonzero)^2 + norm(diag(s))^2  - norm(low_vals_in_nonzero) ^2;
    fro = sqrt(fro);

    % fro2 = norm((A + u*s*v'), 'fro');
end