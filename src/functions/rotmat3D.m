function R = rotmat3D(omega, phi, kappa)
%ROTMAT converts axis (cardinal) and angle to a 3D rotation matrix
%   rotation order: Y-Z-X
%   Input:
%       - omega: rotation around x (scalar) [rad]
%       - phi: rotation around y (scalar) [rad]
%       - kappa: rotation around z (scalar) [rad]
%   Outputs:
%       - R: rotation matrix (3x3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    R = rotmat("x", omega) * rotmat("z", kappa) * rotmat("y", phi);
end
