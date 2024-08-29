function [x_new, v_xy_half] = integrateEquationsOfMotion( ...
    x_old, v_xy_old, E_x_field, B_z_field, time_step, charge_to_mass_ratio ...
)
% INTEGRATEEQUATIONSOFMOTION Update the positions and velocities of the
% particles.
%   This follows section 2-4 of the textbook. The `v_old` is half a time
%   step behind `x_old`.
    % Equation (8).
    v_xy_half = v_xy_old;
    v_xy_half(1) = v_xy_half(1) + charge_to_mass_ratio * E_x_field * (time_step / 2);
    % Equation (9).
    cyclotron_frequency = charge_to_mass_ratio * B_z_field;
    rotation_matrix_x_row = [cos(cyclotron_frequency * time_step); sin(cyclotron_frequency * time_step)];
    rotation_matrix_y_row = [-sin(cyclotron_frequency * time_step); cos(cyclotron_frequency * time_step)];
    v_xy_rotated = v_xy_half;
    v_xy_rotated(1) = rotation_matrix_x_row * v_xy_half;
    v_xy_rotated(2) = rotation_matrix_y_row * v_xy_half;
    % Equation (10).
    v_xy_new = v_xy_rotated;
    v_xy_new(1) = v_xy_rotated(1) + charge_to_mass_ratio * E_x_field * (time_step / 2);

    % Equation (4).
    x_new = x_old + v_xy_new(1) * time_step;
end