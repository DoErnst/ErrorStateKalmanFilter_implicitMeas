function obs = cart2spherical(carts, compForVC)
%CART2SPHERICAL converts given cartesian coordinates to spherical coords
%   according to transformation specified in Velodyne Manual P.55
%   Inputs:
%	    - carts: cartesian coordinates from VLP-16 (nx3)
%	    - compforVC (optional): true, if compensation for vertical correction
%           should be applied
%   Outputs:
%	    - obs: spherical observations as distance, elevation and azimuth (nx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 2
        compForVC = true;
    end
    cart = carts; % avoid overwriting params for function
    if compForVC
        %% compensate for VC
        % if cartesian coordinates include the vertical correction
        % (that they should)
        
        % compute temporary omega values
        tmp_om = asin(cart(:,3) ./ vecnorm(cart,2,2));
        cart(:,3) = applyVC(cart(:,3), tmp_om, false);
    end
    
    obs = nan(size(cart));
    obs(:,1) = vecnorm(cart,2,2);
    obs(:,2) = asin(cart(:,3) ./ obs(:,1));
    obs(:,3) = atan2(cart(:,2), cart(:,1));
    
end

