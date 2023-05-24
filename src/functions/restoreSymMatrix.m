function mat = restoreSymMatrix(vec)
%RESTORESYMMATRIX converts a flattened symmetric matrix back to a matrix
%   Inputs:
%       - vec: flattened symmetric matrix (nx6/21/171)
%   Outputs:
%       - mat: restored matrices (3/6/18x3/6/18xn)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if size(vec, 2) == 6
        n = 3;
    elseif size(vec, 2) == 21
        n = 6;
    elseif size(vec, 2) == 171
        n = 18;
    end
    mat = nan(n, n,size(vec, 1));
    for i = 1 : size(vec, 1)
        tmp = zeros(n);
        for j = 1 : n
            tmp(1:j,j) = vec(i, j*(j+1)/2 - (j-1) :j*(j+1)/2);
        end
        mat(:,:,i) = tmp + tmp' - diag(diag(tmp));
    end
end

