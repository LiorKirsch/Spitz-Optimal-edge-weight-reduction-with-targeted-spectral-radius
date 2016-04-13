function [ y, lambda, it_num ] = power_method_k_sparse_low(A, U, S, V, k, y, it_max, tol ,is_sym,n)
% Find the k largest right eigen values and eigen vectors using the power method.
%
% Input: Z = A + U*S*V' (A is sparse and u,v are basis vector)
%      A - the sparse matrix  (m x n)
%      U - the left low rank basis matrix  (m x k)
%      S - the singular value matrix       (k x k)
%      V - the right low rank basis matrix (n x k)
%      k - the number of eigen values to find
%      y - init vector for the right eigen vector
%      it_max - maximum number of iteration for the power method
%      tol - tolerance for the power method
%      is_sym - (default true) is the matrix Z symetric
%      n - a parmeter that is used to point to the place where needs 
%       to be cut in case the matrix is not symetric  
%

  debug = 0;
  min_iter = 10;


  if ( debug )
    fprintf ( 1, '\n' );
    fprintf ( 1, '     IT      Lambda          Delta-Lambda    Delta-Y\n' );
    fprintf ( 1, '\n' );
  end

  assert( size(y,2) == k, 'y0 should be of size n X k') ; 
  lambda = 0;
  it_num = 0;
  val_dif = 0.0;
  tol_adjsted = tol * norm(A,1);

  ay = A*y;
  if is_sym
     ay = ay + U* S*(V' *y);
  else
     ay(n+1:end,:) = ay(n+1:end,:) + V*S*(U'*y(1:n,:));
     ay(1:n,:) = ay(1:n,:) +  U*S*(V'*y(n+1:end,:));
  end
  

  if k ==1
      lambda = norm ( ay );
      y = ay / lambda;
  else
      [Q,R] = qr(ay,0);
      y = Q;
      lambda = diag(diag(R));
  end
  


  if ( debug )
    fprintf('  %5d  %14e  %14e\n', it_num, lambda, val_dif );
  end

  for it_num = 1 : it_max

      lambda_old = lambda;

      ay = A * y;
      if is_sym
         ay = ay + U* S*(V' *y);
      else
         ay(n+1:end,:) = ay(n+1:end,:) + V*S*(U'*y(1:n,:));
         ay(1:n,:) = ay(1:n,:) +  U*S*(V'*y(n+1:end,:));
      end


      if k ==1
          lambda = norm ( ay );
          y = ay / lambda;
      else
          [Q,R] = qr(ay,0);
          y = Q;
          lambda = diag(diag(R));
      end

      val_dif = norm ( ay - y*lambda,1);


      if ( debug )
         fprintf('  %5d  %14e  %14e \n', it_num, lambda, val_dif);
      end

      if ( val_dif <= tol_adjsted && min_iter <= it_num )
        break
      end 


  end

end