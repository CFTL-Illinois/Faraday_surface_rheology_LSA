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
% Here, A is an 8 by 8 matrix containing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = fhs_general(k, mu, h, q1, q2, Mn)
  % Compute the square of the wave number
  k_squared = k^2;
  
  % ===== SURFACTANT =====
  % Mn is the (n,k)-dependent Marangoni coupling coefficient obtained after
  % eliminating the insoluble-surfactant transport equation:
  %   [[ μ (D^2 + k^2) w_n ]] - Mn * (Dw_n)|_{z=0} = 0
  % Mn = (beta*sigma0*k^2)/(s_n + Ds*k^2)/(1-Gamma0). For clean case, Mn = 0.
  if nargin < 6
    Mn = 0;
  end
  % ===== END SURFACTANT =====

  % Initialize the 8x8 matrix A with zeros
  A = zeros(8,8);

  % Define the z-coordinates for the bottom, interface, and top surfaces
  bottom_surface = -h(1);
  interface_surface = 0;
  top_surface = h(2);

  % Compute exponential terms for the wave number k and q1, q2 at different z-coordinates
  exp_k_z0_plus = exp(  k * bottom_surface);
  exp_k_z1_plus = exp(  k * interface_surface);
  exp_k_z2_plus = exp(  k * top_surface);
  exp_q1_z0_plus = exp( q1 * bottom_surface);
  exp_q1_z1_plus = exp( q1 * interface_surface);
  exp_q2_z1_plus = exp( q2 * interface_surface);
  exp_q2_z2_plus = exp( q2 * top_surface);

  exp_k_z0_minus = exp(  -k * bottom_surface);
  exp_k_z1_minus = exp(  -k * interface_surface);
  exp_k_z2_minus = exp(  -k * top_surface);
  exp_q1_z0_minus = exp( -q1 * bottom_surface);
  exp_q1_z1_minus = exp( -q1 * interface_surface);
  exp_q2_z1_minus = exp( -q2 * interface_surface);
  exp_q2_z2_minus = exp( -q2 * top_surface);

  % Group the exponential terms into vectors for easier manipulation 
  exp_terms_z1_q1 = [exp_k_z1_plus, exp_k_z1_minus, exp_q1_z1_plus, exp_q1_z1_minus];
  exp_terms_z1_q2 = [exp_k_z1_plus, exp_k_z1_minus, exp_q2_z1_plus, exp_q2_z1_minus];
  exp_terms_z0 = [exp_k_z0_plus, exp_k_z0_minus, exp_q1_z0_plus, exp_q1_z0_minus];
  exp_terms_z2 = [exp_k_z2_plus, exp_k_z2_minus, exp_q2_z2_plus, exp_q2_z2_minus];

  % Handle cases where exponential terms become infinite
  if any(isinf(exp_terms_z0))
      exp_terms_z0 = isinf(exp_terms_z0) * 1.0;
  end

  if any(isinf(exp_terms_z2))
      exp_terms_z2 = isinf(exp_terms_z2) * 1.0;
  end

  % Eq1. Kinematic condition. (B1 in KT94)
  % dζ/dt = w @ z=z1
  A(1,1:4) = [ 1,  1,  1,  1] .* exp_terms_z1_q1;

  % Eq2. Continuity of velocity at the interface (B2 in KT94)
  % w1 = w2 @ z=z1
  A(2,1:4) = [ 1,  1,  1,  1] .* exp_terms_z1_q1;
  A(2,5:8) = [-1, -1, -1, -1] .* exp_terms_z1_q2;

  % Eq3. Continuity of velocity derivative at the interface (B3 in KT94)
  % dw1/dz = dw2/dz @ z=z1
  A(3,1:4) = [ k, -k,  q1, -q1] .* exp_terms_z1_q1;
  A(3,5:8) = [-k,  k, -q2,  q2] .* exp_terms_z1_q2;

    % Eq4. Continuity of tangential stresses at the interface (B4 in KT94)
  % μ1(d2/dzz + k2)w1 = μ2(d2/dzz + k2)w2 @ z=z1
  A(4,1:4) =  mu(1) .* [ 2*k_squared, 2*k_squared, (k_squared + q1^2), (k_squared + q1^2)] .* exp_terms_z1_q1;
  A(4,5:8) = -mu(2) .* [ 2*k_squared, 2*k_squared, (k_squared + q2^2), (k_squared + q2^2)] .* exp_terms_z1_q2;
  % ===== SURFACTANT =====
  % Add Marangoni term (after eliminating Γ_n): [[μ(D^2+k^2)w]] - Mn*Dw = 0.
  % Use Dw from layer 1; Eq3 enforces Dw1 = Dw2 at z=0.
  
  %A(4,1:4) = A(4,1:4) + (Mn + k^2 * surface) .* [ k, -k,  q1, -q1] .*
  %exp_terms_z1_q1; %with surface rheology 
  
  A(4,1:4) = A(4,1:4) + Mn.* [ k, -k,  q1, -q1] .* exp_terms_z1_q1;   %without surface rheology
  % ===== END SURFACTANT =====

  % Eq5. No flow condition, bottom surface (B5 in KT94)
  % w1 = 0 @ z=z0
  A(5,1:4) = [ 1,  1,  1,  1] .* exp_terms_z0;

  % Eq6. No slip condition, bottom surface (B6 in KT94)
  % dw1/dz = 0 @ z=z0
  A(6,1:4) = [k, -k, q1, -q1] .* exp_terms_z0;

  % Eq7. No flow condition, top surface (B9 in KT94)
  % w2 = 0 @ z=z2
  A(7,5:8) = [ 1,  1,  1,  1] .* exp_terms_z2;

  % Eq8. No slip condition, top surface (B10 in KT94)
  % dw2/dz = 0 @ z=z2
  A(8,5:8) = [k, -k, q2, -q2] .* exp_terms_z2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
