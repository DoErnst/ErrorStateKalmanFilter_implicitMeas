function F = buildTransitionMat(dT, R, a_hat)
%BUILDTRANSITIONMAT Sets up the transition matrix for the prediction
%   Ref: Sola, 2017 p. 67 (310)
%   Inputs:
%       - dT: time difference to prior epoch [s]
%       - R: rotation matrix to target frame (3x3)
%       - a_hat: corrected accelerations in (a_x,a_y,a_z) [m/s^2]
%   Outputs:
%       - F: transition matrix (16x16)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    F = eye(15);
    F(1:3, 4:6) = eye(3) * dT; % v -> pos
    F(4:6, 7:9) = -axiator(R * a_hat) .* dT; % a -> v
    F(4:6, 10:12) = -R .* dT; % a_b -> v
    F(7:9, 13:15) = -R .* dT; % omega_b -> d_theta

end

