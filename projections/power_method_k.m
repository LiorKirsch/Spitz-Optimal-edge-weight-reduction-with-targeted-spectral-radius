function [ y, lambda, it_num ] = power_method_k(A, k, y, it_max, tol ,stop_crit)

% Find the k largest right eigen values and eigen vectors using the power method.
%
% Input: Z = A + U*S*V' (A is sparse and u,v are basis vector)
%      S - the matrix  (m x n)
%      k - the number of eigen values to find
%      y - init vector for the right eigen vector
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
  
  assert( size(y,2) == k, 'y0 should be of size n X k') ;
  lambda = 0;
  ay = A * y;
  val_dif = inf;
  tol_adjsted = tol * norm(A,1);

  if ( debug )
    fprintf ( 1, '  %5d  %14e  %14e\n', it_num, lambda, val_dif );
  end

  for it_num = 1 : it_max

    lambda_old = lambda;
    ay = A * y ;
    if k==1
      lambda = norm ( ay );
      y = ay / lambda;
    else
      [Q,R] = qr(ay,0);
      y = Q;
      lambda = diag(diag(R));
    end
  
  
    if k ==1 && strcmp(stop_crit, 'norm')
        val_dif = abs(lambda - lambda_old);
    else
        val_dif = norm ( ay - y*lambda,1);
    end
    

    if ( debug )
      fprintf('  %5d  %14e  %14e \n', it_num, lambda, val_dif);
    end

    if ( val_dif <= tol_adjsted && min_iter <= it_num )
      break
    end 

  end

end