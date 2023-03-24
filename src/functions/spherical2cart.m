function cart = spherical2cart(obs, appVC, sc_line)
% SPHERICAL2CART converts given observations into cartesian coords
%   according to transformation specified in Velodyne Manual P.55
%   Inputs:
%   - obs: spherical observations distances, elevation, azimuth (nx3)
%	- appVC: true, if vertical correction should be applied (optional)
%	- sc_line: ring number for points (only for high performance necessary)
%           (optional)
%   Outputs:
%	- cart: cartesian coordinates (nx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License


    if nargin < 2
        appVC = true;
    end

    c2 = cos(obs(:,2,:));
    
    cart = nan(size(obs));
    cart(:,1,:) = obs(:,1,:) .* c2 .* cos(obs(:,3,:));
    cart(:,2,:) = obs(:,1,:) .* c2 .* sin(obs(:,3,:));
    cart(:,3,:) = obs(:,1,:) .* sin(obs(:,2,:));
    
    if appVC
        if nargin < 3
            cart(:,3,:) = applyVC(cart(:,3,:), obs(:,2,:), true);
        else
            cart(:,3,:) = applyVC(cart(:,3,:), obs(:,2,:), true, sc_line);
        end
    end
    
end

