function S = axiator(vecs)
%AXIATOR creates the axiators of the given vectors (works for n-dim)
%   creates the axiator matrices for given vectors
%   Inputs:
%       - vecs: vectors are assumed to be column vectors (3xn)
%   Outputs:
%       - S: axiator matrices (3x3xn)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    n = size(vecs,2);
    vecs = permute(vecs, [1,3,2]);

    S = [zeros(1,1,n), -vecs(3,1,:), vecs(2,1,:);
         vecs(3,1,:), zeros(1,1,n), -vecs(1,1,:);
         -vecs(2,1,:), vecs(1,1,:), zeros(1,1,n)];
end

