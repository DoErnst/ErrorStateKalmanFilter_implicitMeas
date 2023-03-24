function [xd, Sigmaxx, H] = planeEst(coords)
%PLANEEST estimates a plane through a point cloud w/ noise
%   uses eigenvalue decomposition to estimate an optimal plane (least squares)
%   through a point cloud
%   Ref: Drixler: Analyse der Form und Lage von Objekten im Raum, 1993, p.45f
%   Inputs:
%       - coords: cartesian coordinates (nx3)
%   Outputs:
%       - xd: plane params (nx,ny,nz,d; 4x1)
%       - Sigmaxx: VCM of plane params (4x4)
%       - H: trafo matrix to plane-frame (4x4)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    xd = zeros(4,1);
    % params
    A = coords - mean(coords, 1); % (5.3.27)
    [H, Lambda] = eig(A' * A); % (5.3.28)
    xd(1:3) = H(:,1);
    e = ones(size(A,1),1);
    xd(4) = -(e' * e)^-1 * e' * coords * H(:,1); % (5.3.29)
    sigma = sqrt(Lambda(1,1) / (size(coords,1) - 3)); % (5.3.30)
    % VCM
    Qnn = 1 / Lambda(2,2) * H(:,2) * H(:,2)' ...
        + 1 / Lambda(3,3) * H(:,3) * H(:,3)'; % (5.3.31)
    qdd = 1/ size(A,1)^2 * e' * coords * Qnn * coords' * e;
    Qxx = [Qnn, zeros(3,1); zeros(1,3), qdd]; % (5.3.32)
    Sigmaxx = sigma^2 * Qxx;
end