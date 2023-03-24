function plotPcl(pcl,fmtStr)
%PLOTPCL plots point clouds into existing figures
%   Inputs:
%       - pcl: cartesian point cloud (nx3)
%       - fmtStr: formatting string for plot function (string)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 2
        fmtStr = 'b.';
    end

    plot3(pcl(:,1),pcl(:,2),pcl(:,3),fmtStr);
    box on
    grid on
    axis equal

end