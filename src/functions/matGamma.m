function matG = matGamma(vecs)
%MATGAMMA sets up the gamma matrix for the LiDAR simulation
%   Inputs:
%       - vecs: vectors of rays (6xn)
%   Outputs:
%       - matG: gamma matrices (6x6xn)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    n = size(vecs,2);
    S = axiator(vecs(4:6,:));
    vecs = permute(vecs, [1,3,2]);

    matG = [-S, -vecs(1:3,1,:);
            permute(vecs(1:3,1,:),[2,1,3]), zeros(1,1,n)];
end

