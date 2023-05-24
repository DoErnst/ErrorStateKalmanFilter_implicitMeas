function cart = spherical2cart2(obs, altern)
% SPHERICAL2CART converts given observations into cartesian coords
%   conversion for laser tracker/ total station observations (non LiDAR)
%   Inputs:
%       - obs: spherical observations distances, elevation, azimuth (nx3)
%       - altern: swap between definitions for azimuth
%   Outputs:
%       - cart: cartesian coordinates (nx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 2
        altern = false;
    end

    cart = nan(size(obs));
    if altern
        cart(:,1) = obs(:,1) .* sin(obs(:,2)) .* sin(obs(:,3));
        cart(:,2) = obs(:,1) .* sin(obs(:,2)) .* cos(obs(:,3));
    else
        cart(:,1) = obs(:,1) .* sin(obs(:,2)) .* cos(obs(:,3));
        cart(:,2) = obs(:,1) .* sin(obs(:,2)) .* sin(obs(:,3));
    end
    cart(:,3) = obs(:,1) .* cos(obs(:,2));

    
end
    
    