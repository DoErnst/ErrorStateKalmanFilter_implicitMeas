function logIntRegion = pointsInRegion2(pcl, intRegions)
%POINTSINREGION2 checks if points are in given volumes of interest
%   Inputs:
%       - pcl: cartesian coordinates of points (nx3)
%       - intRegion: volume of interest around one panel (6xm)
%   Outputs:
%       - logIntRegion: true if in volume (nxm)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    n_r = size(pcl, 2);

    xcoords = reshape(pcl(1,:,:),[n_r, 16])';
    ycoords = reshape(pcl(2,:,:),[n_r, 16])';
    zcoords = reshape(pcl(3,:,:),[n_r, 16])';

    logMinX = xcoords > intRegions(:,1)';
    logMaxX = xcoords < intRegions(:,2)';
    logMinY = ycoords > intRegions(:,3)';
    logMaxY = ycoords < intRegions(:,4)';
    logMinZ = zcoords > intRegions(:,5)';
    logMaxZ = zcoords < intRegions(:,6)';
    logIntRegion = logMinX & logMaxX & logMinY & logMaxY & logMinZ & logMaxZ;
end

