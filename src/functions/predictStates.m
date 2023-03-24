function [x_k, Q_x_k] = predictStates(x_k_, Q_x_k_, dT, a, omega, s, g_vec)
%PREDICTSTATES Prediction step of the filter based on IMU measurements
%   uses accelerometer and gyroscope measurements to predict states
%   Ref: Sola: Quaternion kinematics for the error-state Kalman filter, 2017 p. 58 (259)
%   Input:
%       - x_k_: states of prior epoch (16x1)
%       - Q_x_k_: VCM of states of prior epoch (16x16)
%       - dT: time difference to prior epoch [s]
%       - a: accelerations in (a_x,a_y,a_z) [m/s^2]
%       - omega: rotation rates around (r_x,r_y,r_z) (3x1) [rad/s]
%       - s: stochastic information/ noise parameters (struct)
%       - g_vec: gravity vector for area of experiment (3x1) [m/s^2]
%   Output:
%       - x_k: predicted states (16x1)
%       - Q_x_k: VCM of predicted states (16x16)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % copy over states
    x_k = x_k_;

    % rotation matrix
    R = rotmat(convertmat2quat(x_k_(7:10)),'point');
    
    % reduce by biases
    a_hat = a - x_k_(11:13);
    omega_hat = omega - x_k_(14:16);

    % rotate accelerometer measurement in global frame and reduce by
    % gravity
    a_glo = R * a_hat - g_vec;

    % predict nominal states

    x_k(1:3) = x_k_(1:3) + x_k_(4:6) * dT + 0.5 * a_glo * dT^2; % pos
    x_k(4:6) = x_k_(4:6) + a_glo * dT; % vel

    d_theta = omega_hat * dT;
    q_d_theta = quaternion(d_theta', 'rotvec');
    x_k(7:10) = convertquat2mat(convertmat2quat(x_k_(7:10)) * q_d_theta);
    x_k(7:10) = x_k(7:10) ./ norm(x_k(7:10));
    % rest of states stay the same

    F = buildTransitionMat(dT, R, a_hat);

    % perturbations
    if size(s.s_a_n,2) > 1
        V_i = diag(s.s_a_n.^2 .* dT^2); % (261)
        Theta_i = diag(s.s_omega_n.^2 .* dT^2); % (262)
        A_i = diag(s.s_a_w.^2 .* dT); % (263)
        Omega_i = diag(s.s_omega_w.^2 .* dT); % (264)
    else
        V_i = s.s_a_n^2 * dT^2 * eye(3); % (261)
        Theta_i = s.s_omega_n^2 * dT^2 * eye(3); % (262)
        A_i = s.s_a_w^2 * dT * eye(3); % (263)
        Omega_i = s.s_omega_w^2 * dT * eye(3); % (264)
    end
    Q_i = blkdiag(V_i, Theta_i, A_i, Omega_i); % (270)
    % jacobian of perturbations
    % isotropic version
    % F_i = zeros(15, 12); % (270)
    % F_i(4:15, 1:12) = eye(12); % (270)

    B = [zeros(3,6); -R, zeros(3); zeros(3), -eye(3); zeros(6)]; % (460)
    C = [zeros(9,6); eye(6)]; % (460)
    F_i = [B, C];

    % prediction
    % d_x = F * d_x; % (267)
    Q_x_k = F * Q_x_k_ * F' + F_i * Q_i * F_i'; % (268)


end

