function mat = convertquat2mat(q)
%CONVERTQUAT2MAT converts a quaternion to a matrix with the same entries
%   Info: used for verbosity in code
%   Input:
%       - q: quaternion (Matlab obj.)
%   Output:
%       - vec: vector with values for quaternion (4x1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License
    
    [qA, qB, qC, qD] = parts(q);
    mat = [qA; qB; qC; qD];

end

