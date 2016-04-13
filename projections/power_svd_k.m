function [ u, s,v, it_num ] = power_svd_k(A, k, u,v, it_max, tol ,stop_crit)

% Find the k largest singular values and singular vectors using the power method.
%
% Input: 
%      A - the matrix  (m x n)
%      k - the number of eigen values to find
%      u - init vector for the left singular vector
%      v - init vector for the right singular vector
%      it_max - maximum number of iteration for the power method
%      tol - tolerance for the power method
%      stop_crit - can be "norm" or "k_eigs" (default: "k_eigs")
%

  debug = 0;
  min_iter = 10;
  if ~exist('stop_crit', 'var')
      stop_crit = 'k_eigs';
  end

  if ( debug )
    fprintf ( 1, '\n' );
    fprintf ( 1, '     IT      Lambda          Delta-Lambda    Delta-Y\n' );
    fprintf ( 1, '\n' );
  end
  
  assert( size(u,2) == k, 'y0 should be of size n X k') ;
  s = 0;
  it_num = 0;
  val_dif = inf;
  A_t = A';
  norm_A = norm(A,1);
  tol_adjusted = tol*norm_A;

  if ( debug )
    fprintf ( 1, '  %5d  %14e  %14e\n', it_num, s, val_dif );
  end

  for it_num = 1 : it_max

    s_old = s;
     if k==1
      av = A * v;
      s = norm ( av );
      u = av / s;

      au = A_t * u;
      s = norm ( au );
      v = au / s;
    else
      av = A * v;
      [Q,R] = qr(av,0);
      u = Q;
      
      au = A_t * u;
      [Q,R] = qr(au,0);
      v = Q;
      s = diag(diag(R));
     end
  
    if k ==1 && strcmp(stop_crit, 'norm')
        val_dif = abs(s - s_old);
    else
        val_dif = norm(au - v*s,1);
    end
  
    if ( debug )
      fprintf('  %5d  %14e  %14e \n', it_num, s, val_dif);
    end

    if ( val_dif <= tol_adjusted && min_iter <= it_num)
      break
    end 

  end

   v = v * sign(s);
   s = abs(s);
end