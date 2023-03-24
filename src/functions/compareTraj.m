function [devs, mds, mds_limits] = compareTraj(trueTraj, filtTraj, covs)
%COMPARETRAJ compares two trajectory (true and some other)
%   computes the deviations between two time synchronized trajectories
%   computations: true - filt
%   the deviations of the trajectories are given as a deviation vector and
%   the deviations of the orientations are rotation vectors between the two
%   orientations
%   Inputs:
%       - trueTraj: [x,y,z,q0,q1,q2,q3]
%       - filtTraj: [x,y,z,.,.,.,q0,q1,q2,q3,...] (nx16)
%       - covs: variance-covariance matrices for poses (15x15xn)
%   Outputs:
%       - devs: deviations in [x,y,z,omega,phi,kappa]
%       - mds: sqrt of Mahalanobis distance of pose
%       - mds_limits: interval of chi2-distribution for pose with 5%
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    n = size(filtTraj,1);

    % deviations in translations
    tra_diffs = trueTraj(:,1:3) - filtTraj(:,1:3);
    rot_diffs = nan(n,3);
    
    % deviations in orientations
    diff_qs = nan(n,4);
    for i = 1 : n
        diff_q = convertmat2quat(trueTraj(i, 4:7)) * conj(convertmat2quat(filtTraj(i, 7:10)));
        diff_qs(i,:) = convertquat2mat(diff_q);
        rot_diffs(i,:) = diff_q.euler('XYZ','point');
    end
    ix_1 = isnan(tra_diffs(:,1));
    ix_2 = isnan(rot_diffs(:,1));
    ix_ = or(ix_1, ix_2);
    tra_diffs(ix_,:) = [];
    rot_diffs(ix_, :) = [];
    devs = [tra_diffs, rot_diffs];
    
    % standard deviations
    stds = nan(n, 6);
    for i = 1 : n
        tempQ = [covs(1:3,1:3,i), covs(1:3,7:9,i); covs(7:9,1:3,i), covs(7:9,7:9,i)];
        stds(i,:) = sqrt(diag(tempQ));
    end
    
    % Mahalanobis distance
    mds = nan(n,1);
    for i = 1 : n
        tempQ = [covs(1:3,1:3,i), covs(1:3,7:9,i); covs(7:9,1:3,i), covs(7:9,7:9,i)];
        mds(i) = devs(i,:) * tempQ^-1 * devs(i,:)';
        mds(i) = sqrt(mds(i));
    end
    mds = mds(~isnan(mds));
    mds_limits = sqrt([chi2inv(0.025, size(devs,2)), chi2inv(0.975, size(devs,2))]); % Bar-Shalom et al., 2004, p. 234
    logInLimits = and(mds > mds_limits(1), mds < mds_limits(2));
    percInLimits = sum(logInLimits) / length(mds);
    disp(['[I] ', num2str(percInLimits * 100, '%.2f'), '% of epochs in interval'])
end

