function R = rotmat(axis, angle)
%ROTMAT converts axis (cardinal) and angle to a 3D rotation matrix
%   Input:
%       - axis: either "x", "y" or "z" (char)
%       - angle: angle for rotation (scalar) [rad]
%   Outputs:
%       - R: rotation matrix (3x3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License
    
    c = cos(angle);
    s = sin(angle);

    if axis == "x" 
        R = [1, 0, 0;
             0, c, -s;
             0, s, c];
    end
    if axis == "y"
        R = [c, 0, s;
             0, 1, 0;
             -s, 0, c];
    end
    if axis == "z"
        R = [c, -s, 0;
             s, c, 0;
             0, 0, 1];
    end
end