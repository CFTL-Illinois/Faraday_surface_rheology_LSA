%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble the full hydrodynamic system matrix
% Taken from the appendix B of Kumar & Tuckerman (1994).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dn] = matrix_An(sigma, alpha, n, omega, k, nu, mu, rho, g, gamma, h, beta, Gamma0, Ds, surface)

  fexp = sigma + 1i * (alpha + n * omega);
  sigma0 = gamma;

  if nargin < 12, beta = 0; end
  if nargin < 13, Gamma0 = 0; end
  if nargin < 14, Ds = 0; end
  if nargin < 15, surface = 0; end

  Mn_M = 0;
  if beta ~= 0 && Gamma0 ~= 0
      Mn_M = (beta * Gamma0 * sigma0 * k^2) / ((fexp + Ds * k^2) * (1 - Gamma0));
  end

  Mn_S = k^2 * surface;
  Mn = Mn_M + Mn_S;

  q1n = sqrt(k.^2 + fexp / nu(1));
  q2n = sqrt(k.^2 + fexp / nu(2));

  interface_surface = 0;
  exp_k_z1_plus = exp(k * interface_surface);
  exp_k_z1_minus = exp(-k * interface_surface);

  eps = 1e-12;

  if fexp ~= 0
    lhs = fhs_general(k, mu, h, q1n, q2n, Mn);

    if (rcond(lhs) > eps)
      rhs = zeros(8, 1);
      rhs(1) = fexp;

      D1 = diag(1 ./ sqrt(sum(abs(lhs), 2)));
      D2 = diag(1 ./ sqrt(sum(abs(lhs'), 2)));
      A_scaled = D1 * lhs * D2;
      b_scaled = D1 * rhs;

      x_scaled = A_scaled \ b_scaled;
      coefficients = D2 * x_scaled;

      exp_q1n_z1_plus = exp( q1n * interface_surface);
      exp_q2n_z1_plus = exp( q2n * interface_surface);
      exp_q1n_z1_minus = exp(-q1n * interface_surface);
      exp_q2n_z1_minus = exp(-q2n * interface_surface);
      exp_terms_z1_q1 = [exp_k_z1_plus, exp_k_z1_minus, exp_q1n_z1_plus, exp_q1n_z1_minus];
      exp_terms_z1_q2 = [exp_k_z1_plus, exp_k_z1_minus, exp_q2n_z1_plus, exp_q2n_z1_minus];

      dz1 = ([k, -k, q1n, -q1n] .* exp_terms_z1_q1) * coefficients(1:4);
      dz2 = ([k, -k, q2n, -q2n] .* exp_terms_z1_q2) * coefficients(5:8);
      dzzz1 = ([k^3, -k^3, q1n^3, -q1n^3] .* exp_terms_z1_q1) * coefficients(1:4);
      dzzz2 = ([k^3, -k^3, q2n^3, -q2n^3] .* exp_terms_z1_q2) * coefficients(5:8);

    else
      lhs = fhs_general_deep(k, mu, q1n, q2n, Mn);

      rhs = zeros(4, 1);
      rhs(1) = fexp;

      D1 = diag(1 ./ sqrt(sum(abs(lhs), 2)));
      D2 = diag(1 ./ sqrt(sum(abs(lhs'), 2)));
      A_scaled = D1 * lhs * D2;
      b_scaled = D1 * rhs;

      x_scaled = A_scaled \ b_scaled;
      coefficients = D2 * x_scaled;

      dz1 = [ k,  q1n] * coefficients(1:2);
      dz2 = [-k, -q2n] * coefficients(3:4);
      dzzz1 = [ k^3, q1n^3 ] * coefficients(1:2);
      dzzz2 = [-k^3, -q2n^3] * coefficients(3:4);
    end

  else
    lhs = fhs_special(k, mu, h, Mn);

    if (rcond(lhs) > eps)
      rhs = zeros(8, 1);
      rhs(1) = fexp;

      D1 = diag(1 ./ sqrt(sum(abs(lhs), 2)));
      D2 = diag(1 ./ sqrt(sum(abs(lhs'), 2)));
      A_scaled = D1 * lhs * D2;
      b_scaled = D1 * rhs;

      x_scaled = A_scaled \ b_scaled;
      coefficients = D2 * x_scaled;

      exp_terms_z1 = [exp_k_z1_plus, exp_k_z1_minus, exp_k_z1_plus, exp_k_z1_minus];
      dz1 = ([k, -k, (1 + k * interface_surface), (1 - k * interface_surface)] .* exp_terms_z1) * coefficients(1:4);
      dz2 = ([k, -k, (1 + k * interface_surface), (1 - k * interface_surface)] .* exp_terms_z1) * coefficients(5:8);
      dzzz1 = ([k^3, -k^3, 3 * k^2 + k^3 * interface_surface, 3 * k^2 - k^3 * interface_surface] .* exp_terms_z1) * coefficients(1:4);
      dzzz2 = ([k^3, -k^3, 3 * k^2 + k^3 * interface_surface, 3 * k^2 - k^3 * interface_surface] .* exp_terms_z1) * coefficients(5:8);

    else
      lhs = fhs_special_deep(k, mu, Mn);

      rhs = zeros(4, 1);
      rhs(1) = fexp;

      D1 = diag(1 ./ sqrt(sum(abs(lhs), 2)));
      D2 = diag(1 ./ sqrt(sum(abs(lhs'), 2)));
      A_scaled = D1 * lhs * D2;
      b_scaled = D1 * rhs;

      x_scaled = A_scaled \ b_scaled;
      coefficients = D2 * x_scaled;

      dz1 = [ k, (1 + k * interface_surface)] * coefficients(1:2);
      dz2 = [-k, (1 - k * interface_surface)] * coefficients(3:4);
      dzzz1 = [ k^3, 3 * k^2 + k^3 * interface_surface] * coefficients(1:2);
      dzzz2 = [-k^3, 3 * k^2 - k^3 * interface_surface] * coefficients(3:4);
    end
  end

  Dn = ((rho(2) * fexp + 3 * mu(2) * k^2) * dz2 - mu(2) * dzzz2 - ...
        (rho(1) * fexp + 3 * mu(1) * k^2) * dz1 + mu(1) * dzzz1 + ...
        ((rho(2) - rho(1)) * g - gamma * k^2) * k^2);
end
