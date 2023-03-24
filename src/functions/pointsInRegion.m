function logIntRegion = pointsInRegion(pcl, intRegion)
%POINTSINREGION checks if points are in given volume of interest
%   Inputs:
%       - pcl: cartesian coordinates of points (nx3)
%       - intRegion: volume of interest around one panel (6x1)
%   Outputs:
%       - logIntRegion: true if in volume (nx1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License 


    logMinX = pcl(:,1) > intRegion(1);
    logMaxX = pcl(:,1) < intRegion(2);
    logMinY = pcl(:,2) > intRegion(3);
    logMaxY = pcl(:,2) < intRegion(4);
    logMinZ = pcl(:,3) > intRegion(5);
    logMaxZ = pcl(:,3) < intRegion(6);
    logIntRegion = logMinX & logMaxX & logMinY & logMaxY & logMinZ & logMaxZ;
end

