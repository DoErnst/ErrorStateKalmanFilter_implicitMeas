function Q_l_k_ = trafoWeightsOfObs(Q_l_k, coords, caliLS2IMU, VCM_LS_IMU)
%TRAFOWEIGHTSOFOBS computes variance propagation for LiDAR observation weights
%   computes the VCM for VCM for the LiDAR observation weights by variance
%   propgation including original (spherical) stddevs and VCM of trafo between
%   LiDAR and IMU
%   Inputs:
%       - Q_l_k: VCM of cartesian points in LiDAR-frame (spherical) (3*nx3*n)
%       - coords: cartesian point coordinates (nx3) [m]
%       - caliLS2IMU: pose of LiDAR in IMU-frame (6x1,tx,ty,tz,omega,phi,kappa)
%       - VCM_LS_IMU: VCM of trafo LiDAR -> IMU (6x6)
%   Outputs:
%       - Q_l_k_: VCM of cartesian coordinates in IMU-frame (3*nx3*n)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    R = quaternion(caliLS2IMU(4:6), 'euler', 'XYZ', 'point').rotmat('point');

    p = size(coords, 1);

    Q_l_k_ = zeros(size(Q_l_k));
    F_rot = buildF_rot(coords, caliLS2IMU);

    for i = 1 : p
        ix = (i-1) * 3;
        F = [R, eye(3), F_rot(ix+1 : ix+3, ix+1 : ix+3)];
        Q_l_k_(ix+1:ix+3,ix+1:ix+3) = F * blkdiag(Q_l_k(ix+1 : ix+3, ix+1 : ix+3), VCM_LS_IMU) * F';
    end

end

