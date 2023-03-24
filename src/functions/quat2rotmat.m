function R = quat2rotmat(q)
%QUAT2ROTMAT converts a quaternion (given as a vector) to a rotation matrix
%   Info: alternative to quaternion(q).rotmat('point') for derivatives
%   Input:
%       - q: quaternion (4x1)
%   Outputs:
%       - R: rotation matrix (3x3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    R = sum(q.^2)^-1 * [
        q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, 2 * (q(2) * q(3) - q(1) * q(4)), 2 * (q(2) * q(4) + q(1) * q(3));
        2 * (q(2) * q(3) + q(1) * q(4)), q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2, 2 * (q(3) * q(4) - q(1) * q(2));
        2 * (q(2) * q(4) - q(1) * q(3)), 2 * (q(3) * q(4) + q(1) * q(2)), q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2
    ];

end

