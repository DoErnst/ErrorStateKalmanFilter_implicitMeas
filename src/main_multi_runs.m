%% Simulations -> for multiple runs
% script for the setup of multiple runs defined in settings file
% results are saved for each trajectory in MAT-files for later analysis
%
% Copyright (c) 2023 Dominik Ernst under MIT License

addpath('functions')
importSimulationInputs; % load settings
disp(['[I] Settings loaded: ', rs.prefix]);
numEpochs = tr.trajLen * fp.IMU_Hz;
iRun = 1;
for t = 1 : rs.numTraj
    %% generate trajectory and sensor data
    if rs.fixRNG
        seed = t;
    else
        t1 = tic;
        seed = mod(t1,2^32);
    end
    rng(seed);
    traj = generateRandomTraj(tr.trajLen, [10,8,6]);
    [IMU_obs, LS_obs, traj_IMU] = generateSensorDataFromTraj(traj, tr, fp);
    fprintf('%s\n',['[I] Trajectory: ', num2str(t, '%i'), ', RNG seed: ', num2str(seed,'%i')]);
    % preallocate matrices
    mat_devs = nan(numEpochs+1, 6, rs.numMCruns); % [m/rad, 6, num MC runs]
    mat_mds = nan(numEpochs+1, rs.numMCruns); % [-, num MC runs]
    times = nan(rs.numMCruns, 1);
    for r = 1 : rs.numMCruns
        %% MC runs
        fprintf('%s',['[I] Run of Trajectory: ', num2str(r, '%i')]);
        tic
        LS_obs_n = addNoiseToPcl(LS_obs, st, fp);
        IMU_obs_n = addNoiseToIMU(IMU_obs, traj_IMU, st, fp);
        % initialize filter
        %   compute accelerometer bias
        R_acc = quaternion(traj_IMU(1,5:8)).rotmat('point');
        acc_b = mean(R_acc * IMU_obs_n(1:300, 2:4)',2) - fp.g_vec;
        %   setup initial state vector
        x0 = [traj(1,2:4), zeros(1,3), traj(1,5:8), acc_b', mean(IMU_obs_n(1:300, 5:7))]';
        Qx0 = diag([0.0001^2 * ones(1,3), 0.001^2 * ones(1,3), deg2rad(0.01)^2 * ones(1,3), ...
               st.s_a_w.^2, st.s_omega_w.^2]);
        % run filter
        [states, statesVCM, times(r)] = ESKF_iM_LiDAR_IMU(LS_obs_n,IMU_obs_n,x0, Qx0, st, fp, pm, tr);
        time_end = toc;
        fprintf('%s\n',[' -> Run took ', num2str(time_end, '%.2f'), ' s.']);
        % compute deviations between true and filter trajectory
        [mat_devs(:,:,r), mat_mds(:,r), ~] = compareTraj(traj_IMU(:,2:end),states, statesVCM);
        if any(any(abs(mat_devs(:,:,r)) > 1))
            disp(['[W] Investigate seed: ', num2str(seed,'%i')]);
        end
    end
    %% save results for further processing
    savePath = '../data/output/';
    nFiles = length(dir(savePath)) - sum([dir(savePath).isdir]) - 1;
    % filename: prefix for identification, number counting up and random 
    %   part to avoid collisions while saving
    fname = [rs.prefix, '_', num2str(nFiles+1,'%04i'),'_',num2str(randperm(999,1),'%03i'),'.mat'];
    save(fullfile(savePath, fname), 'mat_devs', 'mat_mds','times');
    clear mat_devs mat_mds
end
