function [IMU_obs, LS_obs, traj_IMU] = generateSensorDataFromTraj(traj, tr, fp)
%GENERATESENSORDATAFROMTRAJ generates IMU and LiDAR data for given traf.
%   uses a given traj and environment to generate IMU and LiDAR
%   observations
%       Warning: uses a pre-compiled mex-function to faster simulation of
%           point clouds, may need to be re-compiled, alternatively
%           the normal function simulateLiDARdata (slow) or
%           simulateLiDARdata2 (fastest uncompiled) can be used
%   Inputs:
%       - traj: trajectory (t,x,y,z,q0,q1,q2,q3)
%       - tr: info about the planes in env. (struct)
%   Outputs:
%       - IMU_obs: simulated IMU obs. (t,a_x,a_y,a_z,r_x,r_y,r_z)
%       - LS_obs: simulated LiDAR observations in the s-frame
%           (t,x,y,z,sc,obj_id)
%       - traj_IMU: traj interpolated to times of IMU obs.
%           (t,x,y,z,q0,q1,q2,q3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License
    
    tic
    
    [IMU_obs, traj_IMU] = simulateIMUdata(traj, 1/fp.IMU_Hz, eye(4));
    LS_obs = simulateLiDARdata( ...
        traj,eye(4),eye(4),tr.H_r_p,tr.intPP, tr.intRegions, tr.refPts_r, fp.LiDAR_Hz ...
    );
    time_end = toc;
    disp(['[I] Sensor data generated in ', num2str(time_end, '%.2f'), ' s.']);
end

