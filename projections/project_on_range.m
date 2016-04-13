function W_bound_above_and_below =project_on_range(W, A_max)
% returns a matrix (or vector) closest to W and bounded by A_max above and zero below
%
    if isreal(W)
        W_bound_above = min(W, A_max);
        W_bound_above_and_below = max(W_bound_above, zeros(size(A_max)));
    else
%         W_bound_above = min(W, A_max);

        reapart = max(real(W), zeros(size(A_max)));
        imgpart = max(imag(W), zeros(size(A_max)));
        W_bound_above_and_below = complex(reapart,imgpart);
        
    end
end