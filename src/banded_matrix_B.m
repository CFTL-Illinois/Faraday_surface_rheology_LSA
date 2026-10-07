%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function creates B, a banded matrix whose structure depends on alpha
function [B] = banded_matrix_B(N, is_harmonic)
  
  % Assemble Matrix B based on the number of Fourier modes (N) and harmonic flag  
  M = 2*(N+1);

  if is_harmonic == true
    % In the harmonic case, we have
    B = spdiags(ones(M+2,1),  2, M, M);
    B = spdiags(ones(M-2,1), -2, B);
    B(1,3) = 2;
    B(2,4) = 0;
  else
    % and in the subharmonic case, we have
    B = spdiags(ones(M+2,1),  2, M, M);
    B = spdiags(ones(M-2,1), -2, B);
    B(1,1) = 1;
    B(2,2) = -1;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%