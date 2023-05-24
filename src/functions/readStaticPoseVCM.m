function [pose, VCM] = readStaticPoseVCM(fname)
%READSTATICPOSEVCM reads file with a static pose and its VCM and restores them
%   Inputs:
%       - fname: file name with point cloud + additional info
%   Outputs:
%       - pcl: pose (6x1)
%       - VCM: restored VCM (6x6)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    mat = readmatrix(fname);
    pose = mat(1:6)';
    VCM = restoreSymMatrix(mat(7:end));
end

