function mPi = matPi(vec)
%MATPI sets up the Pi matrix for the LiDAR simulation
%   Inputs:
%       - vecs: vector with plane params (4x1)
%   Outputs:
%       - mPi: Pi matrices (4x4)
%
% Copyright (c) 2023 Dominik Ernst under MIT License
    
    mPi = [eye(3) * vec(4), -vec(1:3);
           axiator(vec(1:3)), zeros(3,1)];
end

