function LS_obs = simulateLiDARdata(traj, H_s_t, H_l_r, H_r_p, intPP, intRegions, refPts_r, rr)
%SIMULATELIDARDATA simulates LiDAR observations from traj and planes
%   uses an description of planes in the environment and a provided traj to
%   simulate error-free LiDAR observations of a Velodyne Puck (VLP-16)
%       Info: this function is written to be compiled into a mex-function,
%           as a faster uncompiled version simulateLiDARdata2 can be used
%   Frames:
%       - s: sensor frame of LiDAR
%       - t: frame of true traj
%       - l: local frame of system observing t (often =r)
%       - r: global (room) frame of planes
%   Inputs:
%       - traj: (true) trajectory of system (t,x,y,z,q0,q1,q2,q3)
%       - H_s_t: trafo between LiDAR and (true) traj (4x4)
%       - H_l_r: trafo between local l- and r-frame (4x4)
%       - intPP: plane parameters for j planes in env. (4xj)
%       - intRegions: min and max coords of volumes around planes (6xj)
%       - refPts_r: corner points of planes in r-frame (cell(j))
%       - rr: rotation rate of the LiDAR (5, 10 or 20 Hz, scalar)
%   Outputs:
%       - LS_obs: simulated LiDAR observations in the s-frame
%           (t,x,y,z,sc,obj_id)
%
% Copyright (c) 2023 Dominik Ernst under MIT License


    % general
    n_r = size(refPts_r,1);

    % timestep between columns of VLP-16
    ts = 55.296e-6; % [s]
    if rr == 20
        hor_res = 0.4; % [°]
    elseif rr == 5
        hor_res = 0.1; % [°]
    else
        hor_res = 0.2; % [°]
        if rr ~= 10
            disp('[E] Invalid rotation rate, defaulted to 10 Hz');
        end
    end

    % prepare rays for LiDAR obs
    offsets = [11.2,9.7,8.1,6.6,5.1,3.7,2.2,0.7];
    offsets = [offsets, -offsets(end:-1:1)] * 1e-3;
    tmp_pts = [spherical2cart([ones(16,1), deg2rad(-15:2:15)', zeros(16,1)],false)'; ones(1,16)];
    %tmp_pts = [spherical2cart([ones(16,1), deg2rad(-80:10.66:80)', zeros(16,1)],false)'; ones(1,16)];
    tmp_pts(3,:) = tmp_pts(3,:) + offsets;
    L_s = zeros(6,16);
    for i = 1 : 16
        origin = [0,0,offsets(i),1]';
        L_s(:,i) = matPi(origin) * tmp_pts(:,i);
    end
    % trafo matrix for rotation of scan head
    R_s_s = rotmat3D(0,0,deg2rad(hor_res));
    H_s_s_ = [R_s_s * 1, zeros(3); zeros(3), det(R_s_s)* R_s_s^-1];
    % precompute all directions for scan rays
    L_s_all = zeros(6,16,360/hor_res);
    L_s_all(:,:,1) = L_s;
    for i = 2 : size(L_s_all, 3)
        L_s_all(:,:,i) = H_s_s_ * L_s_all(:,:,i-1);
    end
    timesteps = traj(1,1) : ts : traj(end,1);
    LS_obs = nan(length(timesteps) * 16, 6);
    % interpolate poses for timesteps
    poses_TP = interp1(traj(:,1),traj(:,2:8), timesteps,'pchip');
    poses_TP(:, 4:7) = poses_TP(:, 4:7) ./ vecnorm(poses_TP(:,4:7),2,2);
    for i = 1 : length(timesteps)
        H_t_l = [quat2rotmat(poses_TP(i, 4:7)), poses_TP(i,1:3)'; zeros(1,3),1];
        H_s_r = H_l_r * H_t_l * H_s_t;
        H_r_s = H_s_r^-1;
        % trafo ref planes in s-frame
        Omega_s = (H_r_s^-1)' * intPP;
        % rotation of scan head
        %L_s = H_s_s_ * L_s;
        L_s = L_s_all(:,:, mod(i, size(L_s_all, 3))+1);
        % compute intersection to reference panels
        % [dim, refplane, ring]
        X_s = pagemtimes(matGamma(L_s), Omega_s);
        X_s = X_s ./ X_s(4,:,:);
        X_r = pagemtimes(H_s_r, X_s);  % trafo for boundary check
        X_s = X_s(1:3,:,:);
        % plausibility check using signs of coordinates
        pcheck = sign(X_s) .* permute(sign(L_s(1:3,:)),[1,3,2]);
        X_s(pcheck < 0) = nan;
        % coarse check for points
        % check if points hit within boundaries of panel using ROI
        mask = pointsInRegion2(X_r, intRegions);
        mask3 = permute(repmat(mask, 1, 1, 3), [3,2,1]);
        X_s(~mask3) = nan;
        % add information to points (scan line and hit object)
        refs = permute(repmat(1:n_r,[16,1]),[3,2,1]);
        scanls = permute(repmat(1:16,[n_r,1]),[3,1,2]);
        X_s_ = [X_s; refs; scanls];
        pts_s = reshape(X_s_, 5, 16 * n_r);
        pts_s(:, isnan(pts_s(1,:))) = [];
        % process points by distance to s-origin (prioritize nearer points)
        dists = vecnorm(pts_s(1:3,:),2, 1);
        pts_r = H_s_r * [pts_s(1:3,:); ones(1,size(pts_s,2))];
        [~, s_ix] = sort(dists);
        for j = 1 : length(dists)
            ix = s_ix(j);
            if isnan(pts_s(1,ix)) % skip over already deleted points
                continue
            end
            p_r = H_r_p(:,:,pts_s(4,ix)) * [pts_r(1:3,ix); 1];
            % check if points are in boundary of panel/object
            in_boundary = inPoly(p_r(2:3)', refPts_r{pts_s(4,ix)}(:,2:3));
            if in_boundary
                % all other points of same ray/scanline are further away
                ixs = pts_s(5,:) == pts_s(5,ix);
                ixs(ix) = 0; % don't delete the hit itself
                pts_s(:, ixs) = nan;
            else % no hit
                pts_s(:, ix) = nan;
            end
        end
        % remove columns with nans (non hits)
        pts_s(:, isnan(pts_s(1,:))) = [];
        LS_obs((i-1)*16+1 : (i-1)*16+size(pts_s,2), :) = [timesteps(i) * ones(size(pts_s,2),1), pts_s'];
    end
    LS_obs(isnan(LS_obs(:,2)),:) = [];

end

