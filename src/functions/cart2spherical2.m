function obs = cart2spherical2(cart, altern)
%CART2SPHERICAL converts given cartesian coordinates to spherical coords
%   conversion for laser tracker/ total station observations (non LiDAR)
%   Inputs:
%       - cartesian coordinates (nx3)
%       - altern: swap between definitions for azimuth
%   Outputs:
%	    - obs: spherical observations as distance, elevation and azimuth (nx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 2
        altern = false;
    end

    obs = nan(size(cart));
    obs(:,1) = vecnorm(cart,2,2);
    obs(:,2) = acos(cart(:,3) ./ obs(:,1));
    if altern
        obs(:,3) = atan2(cart(:,1), cart(:,2));
    else
        obs(:,3) = atan2(cart(:,2), cart(:,1));
    end
    
end