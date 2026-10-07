%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The coefficients A, and hence the system (3.12) depend on the Floquet exponent
% σ + iɑ via (3.8) and the boundary conditions (see Appendix B) in a complicated
% manner. However, the amplitude a of the external forcing appears linearly. 
% In fact, for fixed Floquet exponent, (3.12) can be considered to be an 
% eigenvalue problem with eigenvalues 'a' and eigenvectors 'v' whose components 
% are the real and imaginary parts of the ζ_n. 
% That is, we write (3.12) as the generalized eigenvalue problem
%
%   A ζ_n = a B ζ_n
%
% In (3.15), A is a diagonal complex matrix (see matrix_An) and B is a banded 
% matrix whose structure depends on ɑ (see banded_matrix_B).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = compute_eigenvalues(sigma_values, alpha, is_harmonic, k, N, omega, nu, mu, rho, g, gamma, h, beta, Gamma0, Ds, surface)
  nkeep = 50;
  num_sigma_values = length(sigma_values);
  An_block = cell(N + 1, 1);
  a = nan(num_sigma_values, 2 + nkeep);

  % ===== SURFACTANT =====
  % Optional surfactant parameters (default to clean):
  if nargin < 13, beta = 0; end
  if nargin < 14, Gamma0 = 0; end
  if nargin < 15, Ds = 0; end
  if nargin < 16, surface = 0; end
  % ===== END SURFACTANT =====

  % Sweep over temporal frequencies
  for ii = 1:num_sigma_values
    sigma = sigma_values(ii);
    [B] = banded_matrix_B(N, is_harmonic);

    % Sweep over Fourier modes frequencies
    for ij = 1:N + 1
      n = ij - 1;
      An = matrix_An(sigma, alpha, n, omega, k, nu, mu, rho, g, gamma, h, beta, Gamma0, Ds, surface); % SURFACTANT
      An_block{ij} = [real(An) -imag(An); imag(An) real(An)];
    end

    % Assemble the block diagonal matrix A using the matrices in the cell array D
    A = blkdiag(An_block{:});

    % Solve the generalized eigenvalue problem A*x = a*RHS*x.
    % Keep the original KT-style sign convention in RHS, and interpret the
    % physical forcing threshold as abs(real(a)).
    RHS = 0.5 * (rho(2) - rho(1)) * k^2 * g * full(B);
    [~, evals_raw] = eig(A, RHS, "vector");

    % Robust filtering. The old 1e-12 absolute imaginary tolerance often
    % creates artificial holes in neutral curves.
    evals = evals_raw;
    evals = evals(isfinite(evals));

    tolEig = 1e-7;
    evals = evals(abs(imag(evals)) < tolEig * max(1, abs(real(evals))));
    evals = abs(real(evals));

    evals = evals(evals > 1e-10);
    evals = evals(evals < 1e5);
    evals = sort(evals);

    nevals = length(evals);
    a(ii,1) = k;
    a(ii,2) = sigma;
    a(ii,3:2+min(nevals,nkeep)) = evals(1:min(nevals,nkeep));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


