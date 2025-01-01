function [h] = functionChannelEstimates(M, beta, array_spacing, single_side_angle_spread, nom_azi_angle)
    %% Generate channel realizations

    num_paths = 10; % Number of scattering paths

    % Generate \alpha^l \sim \mathcal{N_C}(0,1)
    alpha_l = (randn(num_paths, 1) + 1i*randn(num_paths, 1)) / sqrt(2);

    % Generate \theta^l = \theta + U(-\Delta, \Delta)
    theta_l = nom_azi_angle + (2 * single_side_angle_spread * rand(1, num_paths) - single_side_angle_spread);

    % Initialize the channel vector h
    h = zeros(M, 1);

    % Calculate the contributions from each path
    for l = 1:num_paths
        % Generate the array response vector a(theta^l)
        a_theta_l = exp(1i * 2 * pi * array_spacing * (-(M-1)/2 : (M-1)/2).' * cosd(theta_l(l)));

        % Sum the contributions
        h = h + sqrt(db2pow(beta)) * alpha_l(l) * a_theta_l;
    end

    % Normalize the channel vector
    h = h / sqrt(num_paths);
end
