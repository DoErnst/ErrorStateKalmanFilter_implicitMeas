function matG = matGammaOrigin(vecs)
%MATGAMMAORIGIN sets up the gamma matrix for the LiDAR simulation
%   for simplified computation if the ray pass the origin
%   Inputs:
%       - vecs: vectors of rays (6xn)
%   Outputs:
%       - matG: gamma matrices (6x6xn)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    n = size(vecs,2);
    vecs = permute(vecs, [1,3,2]);

    matG = [zeros(3,3,n), -vecs(1:3,1,:);
            pagetranspose(vecs(1:3,1,:)), zeros(1,1,n)];
end

