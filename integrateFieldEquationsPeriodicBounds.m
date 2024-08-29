function [E_x_field] = integrateFieldEquationsPeriodicBounds( ...
    charge_density, delta_x, N, eps_0 ...
)
% INTEGRATEFIELDEQUATIONSPERIODICBOUNDS Calculate the electric field for
% periodic bounds at `x = 0` and `x = L`.
    charge_density_fft = fft(charge_density);
    wave_numbers = (0:N - 1) / (delta_x * N);
    wave_numbers(wave_numbers >= 1 / (2 * delta_x)) = wave_numbers( ...
        wave_numbers >= 1 / (2 * delta_x) ...
    ) - 1 / delta_x;

    % Section 2-5: equation (12) & (13)
    K = wave_numbers * (sin(wave_numbers * delta_x / 2) / (wave_numbers * delta_x / 2));
    potential_fft = charge_density_fft / (eps_0 * K^2);
    % Section 2-5: equaiton (10) & (11)
    kappa = k * (sin(wave_numbers * delta_x) / (wave_numbers * delta_x));
    E_x_field_fft = -1i * kappa * potential_fft;

    E_x_field = ifft(E_x_field_fft);
end