function [expIMU, traj_int, traj_IMU] = simulateIMUdata(traj,delta_t, H_i_t, movmeanFac)
%SIMULATEIMUDATA simulates IMU observations based on a provided traj
%   is based on the provided delta_t for the sampling rate
%   can include an additional transformation to the true trajectory
%   Inputs:
%       - traj: trajectory (t,x,y,z,q0,q1,q2,q3)
%       - delta_t: time differences in [s]
%       - H_i_t: trafo H-matrix between frame of IMU and (true) traj (4x4)
%       - movmeanFac: factor for moving mean smmothing (scalar, default=1 (off))
%   Outputs:
%       - expIMU: simulated IMU obs. (t,a_x,a_y,a_z,r_x,r_y,r_z)
%       - traj_int: interpolated traj of system (t,x,y,z,q0,q1,q2,q3)
%       - traj_IMU: interpolated traj of IMU (t,x,y,z,q0,q1,q2,q3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 4
        movmeanFac = 1;
    end

    % interpolate original trajectory
    timesteps = traj(1,1) : delta_t : traj(end, 1);
    traj_int = [timesteps', interp1(traj(:,1), traj(:,2:8),timesteps,'pchip')];
    traj_int(:, 5:8) = traj_int(:, 5:8) ./ vecnorm(traj_int(:, 5:8),2,2);
    traj_IMU = [traj_int(:,1), zeros(size(traj_int,1),7)];
    % transform from l-frame to i-frame
    rot_delta_t = zeros(size(traj_int, 1),3);
    t_delta_t = zeros(size(traj_int, 1),3);
    H_t_l = [quat2rotmat(traj_int(1, 5:8)), traj_int(1,2:4)'; 0,0,0,1];
    H_l_i_old = (H_t_l * H_i_t)^-1;
    H_i_l = H_l_i_old^-1;
    traj_IMU(1,2:4) = H_i_l(1:3,4);
    traj_IMU(1,5:8) = quaternion(H_i_l(1:3,1:3),'rotmat','point').compact();
    for i = 2 : size(traj_int, 1)
        H_t_l = [quat2rotmat(traj_int(i, 5:8)), traj_int(i,2:4)'; 0,0,0,1];
        H_l_i_new = (H_t_l * H_i_t)^-1;
        H_i_l = H_l_i_new^-1;
        traj_IMU(i,2:4) = H_i_l(1:3,4);
        traj_IMU(i,5:8) = quaternion(H_i_l(1:3,1:3),'rotmat','point').compact();
        H_diff = H_l_i_old * H_l_i_new^-1;
        t_delta_t(i, 1:3) = H_diff(1:3,4) ./ delta_t;
        rot_delta_t(i, 1:3) = quaternion(H_diff(1:3,1:3),'rotmat','point').rotvec() ./ delta_t;
        H_l_i_old = H_l_i_new;
    end 
    v_delta_t = [zeros(1,3); diff(t_delta_t, 1, 1) ./ delta_t];
    IMU_obs_LT = [timesteps', v_delta_t, rot_delta_t];
    expIMU = [IMU_obs_LT(:,1), movmean(IMU_obs_LT(:,2:7),movmeanFac)];

end

