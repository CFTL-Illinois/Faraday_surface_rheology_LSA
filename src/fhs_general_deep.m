%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble the full hydrodynamic system matrix for the general case
% Taken from the appendix B of Kumar & Tuckerman (1994).
%
% In the Floquet expansions (3.1) and (3.2), the nth coefficient w_jn of the
% vertical velocity in the jth fluid layer is given by (3.7):
%
%  w_jn(z)= a_jn exp(kz) + b_jn exp(-kz) + c_jn exp(q_jn z) + d_jn exp(-q_jn z)
%
% The eight coefficients in (3.7) are expressed in terms of the nth coefficient
% of the interface position by means of the eight equations (2.27) and
% (2.19)-(2.25).
%
% Here, A is an 4 by 4 matrix containing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = fhs_general_deep(k, mu, q1, q2, Mn)
  % Compute the square of the wave number
  k_squared = k^2;

  % ===== SURFACTANT =====
  if nargin < 5
    Mn = 0;
  end
  % ===== END SURFACTANT =====

  % Initialize the 4x4 matrix A with zeros
  A = zeros(4,4);

  % Eq1. Kinematic condition. (B1 in KT94)
  % dζ/dt = w @ z=z1
  A(1,:) = [ 1,  1,  0,  0];

  % Eq2. Continuity of velocity at the interface (B2 in KT94)
  % w1 = w2 @ z=z1
  A(2,:) = [ 1,  1, -1, -1];

  % Eq3. Continuity of velocity derivative at the interface (B3 in KT94)
  % dw1/dz = dw2/dz @ z=z1
  A(3,:) = [ k, q1,  k, q2];

  % Eq4. Continuity of tangential stresses at the interface (B4 in KT94)
  % μ1(d2/dzz + k2)w1 = μ2(d2/dzz + k2)w2 @ z=z1
  A(4,1:2) =  mu(1) .* [ 2*k_squared, (k_squared + q1^2)];
  A(4,3:4) = -mu(2) .* [ 2*k_squared, (k_squared + q2^2)];
  % ===== SURFACTANT =====
  % [[μ(D^2+k^2)w]] - Mn*Dw = 0, using Dw from layer 1 at z=0.
  A(4,1:2) = A(4,1:2) + Mn .* [ k,  q1];
  % ===== END SURFACTANT =====

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
