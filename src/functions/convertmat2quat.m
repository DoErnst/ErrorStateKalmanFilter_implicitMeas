function q = convertmat2quat(vec)
%CONVERTMAT2QUAT converts a vector to a quaterion
%   Info: used for verbosity in code
%   Input:
%       - vec: vector with values for quaternion (4x1)
%   Output:
%       - q: quaternion (Matlab obj.)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if size(vec, 2) < 4
        vec = vec';
    end
    q = quaternion(vec);

end

