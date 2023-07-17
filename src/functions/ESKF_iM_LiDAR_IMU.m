function [states,covs,time] = ESKF_iM_LiDAR_IMU(LS_obs, IMU_obs, x0, Qx0, st, fp, pm, tr)
%ESKF_iM_LiDAR_IMU executes a Error State Kalman Filter with implicit m.eq.
%   uses LiDAR and IMU observations for pose estimation
%   IMU is used for prediction and LiDAR pts are assigned to known planes
%   Refs: Sola: Quaternion kinematics for the error-state Kalman filter, 2017 p. 58 (259)
%         Vogel: Iterated Extended Kalman Filter with Implicit Measurement Equation and
%           Nonlinear Constraints for Information-Based Georeferencing, 2018
%   Inputs:
%       - LS_obs: LiDAR pts in s-frame of LiDAR (t,x,y,z)
%       - IMU_obs: accelerations and rotation rates
%           (t,a_x,a_y,a_z,r_x,r_y,r_z)
%       - x0: initial state vector (16x1)
%       - Qx0: initial VCM of states (16x16)
%       - st: stochastic model for sensors (struct)
%       - fp: functional parameters (trafos, dist offset, etc) (struct)
%       - pm: parameters for the filter (struct)
%       - tr: info about the planes in env. (struct)
%   Outputs:
%       - states: state vector for each epoch of IMU_obs (nx16)
%       - covs: VCM of states for each epoch of IMU_obs (16x16xn)
%       - times: run times of filter [s]
%
% Copyright (c) 2023 Dominik Ernst under MIT License
    
    % check inputs -> set defaults otherwise (some with warnings)
    % also acts as documentation for the necessary filter inputs
    [st, fp, pm] = checkFilterInputs(st, fp, pm);

    % Sola, 2017 sec: 7.2
    % nominal state vector -> error states with prefix d_
    % p: position in global frame (3x1) [m]
    % v: velocity wrt global frame (3x1) [m/s]
    % q: orientation wrt global frame (4x1/ 3x1 in error state) [-]
    % a_b: accelerometer bias (local) (3x1) [m/s^2]
    % omega_b: gyroscope bias (local) (3x1) [rad/s]
    
    %% Initialization
    u = 16; % dimension of full state vector
    % start loop for multiple (random) runs
    states = nan(length(IMU_obs),u);
    covs = nan(u-1,u-1,length(IMU_obs));
    % initial state and VCM
    x_k = x0;
    Q_x_k = Qx0;
    x_k_old = x_k;
    % iteration over all IMU observations
    states(1,:) = x_k;
    k = 1; % index IMU obs
    currentTime = 0; % [s]
    pcl_ix = 1;
    fp.H_LS_IMU = buildH_mass(fp.caliLS2MS,1,'XYZ');
    %% Start Iteration
    tic
    while k <= length(IMU_obs)
        
        % time difference
        dT = IMU_obs(k, 1) - currentTime;
    
        % measurements of the epoch
        a = fp.R_corr_IMU * IMU_obs(k, 2:4)'; % [m/s^2] 
        omega = IMU_obs(k, 5:7)'; % [rad/s]
    
        %% Prediction
        [x_k, Q_x_k] = predictStates(x_k, Q_x_k, dT, a, omega, st, fp.g_vec);
        
        %% Update
        if pcl_ix < size(LS_obs,1) % are still LiDAR obs left?
            %% Processing of LiDAR points           
            % find points for current epoch
            % pcl_time as time of mean pcl
            %new_pcl_ix = find(LS_obs(:,1) > IMU_obs(k,1), 1, 'first');
            new_pcl_ix = pcl_ix;
            while LS_obs(new_pcl_ix,1) < IMU_obs(k,1)
                new_pcl_ix = new_pcl_ix + 1;
                if new_pcl_ix > size(LS_obs,1)
                    new_pcl_ix = new_pcl_ix - 2;
                    break
                end
            end
            curr_pcl = LS_obs(pcl_ix:new_pcl_ix,:);
    
            % apply distance offset
            curr_pcl(:,2:4) = spherical2cart( ...
                cart2spherical(curr_pcl(:,2:4)) + [ones(size(curr_pcl,1),1) * fp.r0_LS, zeros(size(curr_pcl,1),2)] ...
            );
            
            %% transform pcl into IMU frame
            % origin of observations in LS, but pose relating to IMU
            % variance propagation for transformation from LS to IMU to
            %   achieve better Qll for observations without introducing more
            %   complicated derivations
            pcl_ = [curr_pcl(:,2:4)';ones(1, size(curr_pcl,1))];
            H_LS2MS = buildH_mass(fp.caliLS2MS,1,'XYZ');
            tpcl = H_LS2MS * pcl_; % -> work with these coordinates
                        
            % transformation of point cloud using current pose
            H_k = [convertmat2quat(x_k(7:10)).rotmat('point'), x_k(1:3); 0,0,0,1];
            rpcl = H_k * tpcl;
        
            %% assignment of points to panels
            % >matrix with observations and corresponding plane parameters
            % asgdPlanes = [ax, ay, az, nd] for each observation (w/ repetitions)
            % l_k = [x,y,z,x,y, ...] in [m]
            [coords_ix, asgdPlanes] = findPtsInROIforFilter(tr, rpcl, pm.REM_OUTLIERS_LS);
            % uncomment next 2 lines if captured plane indices are known
            %coords_ix = transpose(1:size(rpcl,2));
            %asgdPlanes = tr.intPP(:, curr_pcl(:,5))';
            if ~isempty(coords_ix) % only perform update if points are assigned
                % reduce number of points for faster update -> subsample
                % problem: subsampling of the whole cloud might lead to
                %   discarding points of small planes losing important
                %   information
                % experimental solution here: subsample per plane
                %   pm.MIN_PTS_PER_PLANE sets the minimum number of pts
                %   default of 0 means normal subsampling
                if pm.SUBSAMPLING_FACTOR < 1
                    % get larger from subsampled points or minimum pts
                    %   minimum pts gets checked against available pts
                    n_sub = max(ceil(length(coords_ix) * pm.SUBSAMPLING_FACTOR), ...
                                min(pm.MIN_PTS_PER_UPDATE, length(coords_ix)));
                    r_ix = randsample(size(coords_ix,1), n_sub);
                    r_ix = sort(r_ix);
                    coords_ix = coords_ix(r_ix);
                    asgdPlanes = asgdPlanes(r_ix, :);
                end
                coords_k = tpcl(1:3, coords_ix)';
                times_k = curr_pcl(coords_ix, 1);
                p = size(coords_k, 1);
                n = p * 3;
                l_k = reshape(coords_k', n, 1);
                %% stochastic model
                sigma2_LS = diag(st.VCM_LS);
                % extend by uncertainty of dist offset
                sigma2_LS(1) = sigma2_LS(1) + st.std_r0_LS^2;
                % repetition of observation variances (spherical)
                Q_l_s = sparse(diag(repmat(sigma2_LS, p, 1)));
                % propagation spherical -> cartesian (s-frame)
                F_vp = sparse(buildF_s2c(curr_pcl(coords_ix,2:4)));
                Q_l_k = F_vp * Q_l_s * F_vp';
                % propagation s-frame -> b-frame
                Q_l_k = trafoWeightsOfObs(Q_l_k, curr_pcl(coords_ix,2:4), fp.caliLS2MS, st.VCM_LS_IMU);
                Q_l_k = sparse(Q_l_k); % for safety
                % time handling for interpolation
                if pm.NO_LERP
                    tau = ones(p, 1);
                else
                    tau = (times_k - currentTime) ./ (dT);
                    tau(isnan(tau)) = 1;
                end
                %% update iteration
                % standard formulas from Vogel, 2018: (7) - (19)
                x_k_m = x_k;
                l_k_m = l_k;
                old_d_x_k = inf * ones(15,1);
                for i = 1 : 5
                    coords_k_m = reshape(l_k_m', 3, p)';
                    % computation of design (w/Xdx), restriction matrix and
                    % misclosure vector
                    h = compute_h(x_k_m, coords_k_m, asgdPlanes, x_k_old, tau);
                    A = buildA(x_k_m, coords_k_m, asgdPlanes, x_k_old, tau); % (7 -> changed)
                    B = sparse(buildB(x_k_m, coords_k_m, asgdPlanes, x_k_old, tau)); % (8)
                    % computation state difference (x_km - x_k) (rotation with rotvecs)
                    delta_s = zeros(15,1);
                    delta_s(1:6) = x_k(1:6) - x_k_m(1:6);
                    d_theta = convertmat2quat(x_k(7:10)) * conj(convertmat2quat(x_k_m(7:10)));
                    delta_s(7:9) = d_theta.rotvec(); % ori
                    delta_s(10:15) = x_k(11:16) - x_k_m(11:16);
                    O = A * Q_x_k * A'; % (9)
                    S = full(B * Q_l_k * B'); % (10)
                    invQ_r = (O + S)^-1;
                    K = Q_x_k * A' * invQ_r; % (11)
                    r = B * (l_k - l_k_m) + A * delta_s; % (12)
                    d_x_k = -K * (h + r); % compute error state after (13)
                    
                    % inject error state
                    x_k_m(1:6) = x_k(1:6) + d_x_k(1:6);
                    d_theta = quaternion(d_x_k(7:9)','rotvec');
                    x_k_m(7:10) = convertquat2mat(d_theta * convertmat2quat(x_k(7:10)));
                    x_k_m(7:10) = x_k_m(7:10) ./ norm(x_k_m(7:10));
                    x_k_m(11:16) = x_k(11:16) + d_x_k(10:15);
    
                    if norm(d_x_k - old_d_x_k) < 1e-4
                        break
                    end
                    old_d_x_k = d_x_k;
    
                    l_k_m = l_k - (Q_l_k * B' * invQ_r * (h + r)); % or invQ_r -> S^-1 after Dang, 2008 (14+15)
                end
                x_k = x_k_m;
                % update to cov of states
                L = eye(15) - K * A; % (18)
                Q_x_k = L * Q_x_k * L' + K * S * K'; % (19)
                % additional update to cov of states for quat
                G = eye(15); % (287)
                G(7:9,7:9) = eye(3) + 0.5 .* axiator(d_x_k(7:9)); % (287)
                Q_x_k = G * Q_x_k * G'; % (285)
            end 
        end
        %% finish of epoch
        x_k_old = x_k;
        % save results
        states(k,:) = x_k;
        covs(:,:,k) = Q_x_k;
        % increment counters
        k = k + 1;
        currentTime = currentTime + dT;
        pcl_ix = new_pcl_ix + 1;
    end
    time = toc;

end