function plotDevsMds(times,devs,mds,rs)
%PLOTDEVS Plots the deviations of the trajectory 
%   creates a deviation, RMSE and Mahalanobis distance plot
%   Inputs:
%       - times: vector with times of trajectory (nx1) [s]
%       - devs: matrix with deviations (nx6xm) [m|rad]
%       - mds: matrix with Mahalanobis distances (nxm) [-]
%       - rs: settings for multiple runs (struct)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 4
        MCruns = 1;
    else
        MCruns = rs.numMCruns;
    end

    times = reshape(times, 1, length(times));

    %% deviations
    figure
    meanDevs = mean(devs, 3);
    % translations
    subplot(2,1,1)
    plot(times,meanDevs(:,1) .* 1e3,'Color',[1,0,0])
    hold on
    plot(times,meanDevs(:,2) .* 1e3,'Color',[0,1,0])
    plot(times,meanDevs(:,3) .* 1e3,'Color',[0,0,1])
    box on
    grid on
    xlabel('time [s]')
    ylabel('mean deviation [mm]')
    title("Translations")
    xlim([times(1), times(end)])
    ylim([min(min(meanDevs(:,1:3).*1e3)),max(max(meanDevs(:,1:3).*1e3))])
    % orientations
    subplot(2,1,2)
    plot(times,rad2deg(meanDevs(:,4)) .* 1e3,'Color',[1,0,0])
    hold on
    plot(times,rad2deg(meanDevs(:,5)) .* 1e3,'Color',[0,1,0])
    plot(times,rad2deg(meanDevs(:,6)) .* 1e3,'Color',[0,0,1])
    box on
    grid on
    xlabel('time [s]')
    ylabel('mean deviation [mdeg]')
    title("Orientations")
    xlim([times(1), times(end)])
    ylim([min(min(rad2deg(meanDevs(:,4:6).*1e3))),max(max(rad2deg(meanDevs(:,4:6).*1e3)))])
    %% RMSE
    figure
    rmse_tra = sqrt(mean(mean(devs(:,1:3,:).^2,2),3));
    rmse_rot = sqrt(mean(mean(devs(:,4:6,:).^2,2),3));
    rmse_tot = sqrt(mean(mean(devs(:,1:6,:).^2,2),3));
    rmse_tra_full = sqrt(mean(rmse_tra.^2));
    rmse_rot_full = sqrt(mean(rmse_rot.^2));
    rmse_tot_full = sqrt(mean(rmse_tot.^2));
    disp(['[I] RMSE: ', num2str(rmse_tra_full*1e3,'%.2f'),' mm , ',num2str(rad2deg(rmse_rot_full)*1e3,'%.2f'),' mdeg']);
    disp(['[I] RMSE: ', num2str(rmse_tot_full*1e3,'%.2f'),' [mm|rad]']);
    % translations
    subplot(2,1,1)
    plot(times,rmse_tra .* 1e3,'Color',[0,0,1])
    box on
    grid on
    xlabel('time [s]')
    ylabel('RMSE [mm]')
    title("Translations")
    xlim([times(1), times(end)])
    ylim([0,max(max(rmse_tra) .* 1e3)])
    % orientations
    subplot(2,1,2)
    plot(times,rad2deg(rmse_rot) .* 1e3,'Color',[0,0,1])
    hold on
    box on
    grid on
    xlabel('time [s]')
    ylabel('RMSE [mdeg]')
    title("Orientations")
    xlim([times(1), times(end)])
    ylim([0,max(max(rad2deg(rmse_rot) .* 1e3))])
    %% Mahalanobis distance
    figure
    meanMds = mean(mds,2);
    medianMds = median(mds,2);
    mmds = squeeze(mean(reshape(mds,size(mds,1),size(mds,2)/MCruns,size(mds,2)/(size(mds,2)/MCruns)),2));
    minMds = min(mmds, [], 2);
    maxMds = max(mmds, [], 2);
    plot(times, meanMds(:,1),'Color',[0,0,1],'LineWidth',0.1)
    hold on
    if size(mds, 2) > 1
        plot(times, medianMds(:,1),'Color',[0,0,0],'LineWidth',0.1)
        fill([times, times(end:-1:1)],[minMds; maxMds(end:-1:1)],'b','FaceAlpha',0.3,'EdgeColor','none')
    end
    % interval lines
    mds_limits = sqrt(...
        [chi2inv(0.025, size(devs,2)*MCruns), ...
         chi2inv(0.975, size(devs,2)*MCruns)]./MCruns); % Bar-Shalom et al., 2004, p. 234
    line([times(1), times(end)],[mds_limits(1),mds_limits(1)],'Color',[0,1,0])
    line([times(1), times(end)],[mds_limits(2),mds_limits(2)],'Color',[0,1,0])
    hold off
    box on
    grid on
    xlabel('time [s]')
    ylabel('sqrt(d_m) [-]')
    title("Square Root of Mahalanobis distance")
    xlim([times(1), times(end)])
    ylim([max(min(minMds),1e-0),max([maxMds; mds_limits(2)])])
    if max(maxMds) > 10 * mds_limits(2)
        set(gca, 'YScale', 'log')
    end
end

