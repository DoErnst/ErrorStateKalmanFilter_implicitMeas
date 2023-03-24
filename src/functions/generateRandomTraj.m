function traj = generateRandomTraj(len, roomSize)
%GENERATERANDOMTRAJ generates random poses over time in the given room
%   takes a length in seconds and a room size to generate random poses
%   afterwards, the time between poses is determined to be long enough to
%   limit the velocities in translations and rotations to be within
%   specified values
%   Inputs:
%       - len: traj length [s]
%       - roomSize: (x,y,z) starting at (0,0,0) [m]
%   Output:
%       - traj: trajectory (t,x,y,z,q0,q1,q2,q3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % CONSTANTS to limit dynamics of system
    MAX_TRA_VEL = 1; % m/s
    MAX_ROT_VEL = .4; % rad/s
    STATIONARY_TIME = 5; % stationary time at the beginning of the traj

    % generate poses (1 per second = too many, but gets shortened later)
    traj_ = generateRandomPoses(len, roomSize);
    %% add timestamps to poses
    % limit translational and/or rotational speed
    tra_diffs = vecnorm(traj_(1:end-1,1:3) - traj_(2:end,1:3),2,2);
    time_diffs = tra_diffs ./ MAX_TRA_VEL;
    rot_diffs = quaternion(traj_(1:end-1,4:7)) .* quaternion(traj_(2:end,4:7));
    quatMat = rot_diffs.compact;
    diff_rv = 2 .* acos(quatMat(:,1));
    time_diffs = max(time_diffs, diff_rv ./ MAX_ROT_VEL);
    % fill in additional poses at start (repeat to get a stationary period)
    traj_ = [traj_(1,:);traj_(1,:);traj_];
    ts = [0; cumsum([STATIONARY_TIME/2;STATIONARY_TIME/2; time_diffs])];
    traj_ = [ts, traj_]; % full, too long traj
    %% check dynamics of trajectory again, since check above corresponds to
    %   linear interpolation and later we use pchip
    timesteps = traj_(1,1) : 0.01 : traj_(end, 1);
    traj_int = interp1(traj_(:,1), traj_(:,2:8),timesteps,'pchip');
    traj_int(:, 4:7) = traj_int(:, 4:7) ./ vecnorm(traj_int(:, 4:7),2,2);
    % convert quaternions to angular velocities
    rvs = quaternion(traj_int(:,4:7)).angvel(0.01, 'point');
    rvs(1,:) = [0,0,0]; % otherwise first entry is way too large
    % compute velocities
    tra_vels = [0;vecnorm(diff(traj_int(:, 1:3), 1, 1) ./ 0.01,2,2)];
    rot_vels = [0;max(abs(rvs),[],2)];
    ts_new = ts; % make a copy so the original times are unchanged in loop
    for i = 2 : length(ts)
        ixs = and(timesteps > ts(i-1), timesteps < ts(i));
        time_factor = max(max(tra_vels(ixs)) ./ MAX_TRA_VEL, ...
                       max(rot_vels(ixs)) ./ MAX_ROT_VEL);
        if time_factor > 1
            ts_new(i : end) = ts_new(i : end) + (time_factor * (ts(i)-ts(i-1)));
        end
    end
    traj_ = [ts_new, traj_(:,2:end)]; % full, too long traj
    %% get only poses up to time length of traj
    times = sort([traj_(:,1);len]);
    times = times(times <= len);
    traj = [times, interp1(traj_(:,1),traj_(:,2:end), times, 'pchip')];
    % safety: normalize quaternions
    traj(:,5:8) = traj(:,5:8) ./ vecnorm(traj(:,5:8),2,2);
end