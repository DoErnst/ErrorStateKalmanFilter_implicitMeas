function IMU_obs = addNoiseToIMU(IMU_obs_LT, traj_IMU, st, fp)
%ADDNOISETOIMU adds measurement noise to IMU measurements
%   uses simulated error-free IMU observations and adds the specified noise
%   Inputs:
%       - IMU_obs_LT: simulated accelerations and rotation rates
%           (t,a_x,a_y,a_z,r_x,r_y,r_z)
%       - traj_IMU: trajectory of simulated IMU (t,x,y,z,q0,q1,q2,q3)
%       - st: struct with stochastic information (noise parameters)
%   Outputs:
%       - expIMU: IMU observations with added noise 
%           (t,a_x,a_y,a_z,r_x,r_y,r_z)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % time between measurements
    dt = mean(diff(IMU_obs_LT(:,1),1,1));
    noise = randn(size(IMU_obs_LT(:,2:7))) .* [st.s_a_n, st.s_omega_n] ./ sqrt(dt);
    bias_i = randn(size(IMU_obs_LT(:,2:7))) .* [st.s_a_w, st.s_omega_w] .* sqrt(dt);
    IMU_obs = IMU_obs_LT;
    IMU_obs(:,2:7) = IMU_obs(:,2:7) + noise + bias_i;
    % add gravity acc to acc meas
    rots = quaternion(traj_IMU(:,5:8)).conj.rotmat('point');
    gs = repmat(fp.g_vec, 1, 1, size(IMU_obs_LT, 1));
    gs_s = reshape(pagemtimes(rots, gs), 3, size(IMU_obs_LT, 1))';
    IMU_obs(:,2:4) = IMU_obs(:,2:4) + gs_s;
end

