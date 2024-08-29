function [E_x_field] = integrateFieldEquationsElectrodeBounds( ...
    charge_density, voltage_at_0, voltage_at_L, delta_x, N, eps_0 ...
)
% INTEGRATEFIELDEQUATIONSELECTRODEBOUNDS Calculate the electric field given
% the potential at `x = 0` and `x = L`.
%   This section uses equation (3) of appendix D
    % Appendix D, equation (3).
    right_hand_vector = -charge_density * delta_x^2 / eps_0;
    right_hand_vector(1) = right_hand_vector(1) - voltage_at_0;
    right_hand_vector(end) = right_hand_vector(end) - voltage_at_L;

    derivative_matrix = spdiags([-1 2 -1], -1:1, N - 1, N - 1);
    potential = derivative_matrix\right_hand_vector;
    
    % Section 2-5, equation (4).
    E_x_field = zeros(N - 1);
    E_x_field(1) = (voltage_at_0 - potential(2)) / (2 * delta_x);
    E_x_field(N - 1) = (potential(N - 2) - voltage_at_L) / (2 * delta_x);
    E_x_field(2:N - 2) = (potential(1:N - 3) - potential(3:N - 1)) / (2 * delta_x);
end
