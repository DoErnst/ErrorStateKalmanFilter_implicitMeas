function F = buildF_s2c(coords)
%BUILDF_S2C jacobian matrix for the var.prop. spherical to cartesian
%   Inputs:
%       - coords: cartesian coordinates in s-frame (nx3)
%   Outputs:
%       - F: jacobian matrix (n*3xn*3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    obs = cart2spherical(coords);
    di = obs(:,1)';
    el = obs(:,2)';
    az = obs(:,3)';


    tempF = [cos(az) .* cos(el);
            -di .* cos(az) .* sin(el);
            -di .* cos(el) .* sin(az);
            cos(el) .* sin(az);
            -di .* sin(az) .* sin(el);
             di .* cos(az) .* cos(el);
            sin(el);
            di .* cos(el);
            zeros(1, size(coords,1))];
    
    colF = reshape(tempF, 3, size(coords, 1) * 3)';

    F = zeros(size(coords, 1) * 3, size(coords, 1) * 3);
    for i = 1 : size(coords, 1)
        F((i-1)*3+1 : (i-1)*3+3, (i-1)*3+1 : (i-1)*3+3) = colF((i-1)*3+1 : (i-1)*3+3,:);
    end
end

