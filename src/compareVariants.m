%% Compare variants
% script to compare different variants of the filter
% parses through results and picks relevant results by prefix
%
% Copyright (c) 2023 Dominik Ernst under MIT License

addpath('functions')
importSimulationInputs;
PATH = "..\data\output\";
% set variants to compare:
vars = {'base', 'nolerp'};
n = length(vars);
% maximum runs for comparison (inf = minimum # of files from output)
MAX_FILES = 10;
% preallocate
d = dir(PATH);
d = d(~ismember({d.name},{'.','..'}));
dev_cell = cell(n, length(d));
mds_cell = cell(n, length(d));
ffs = zeros(n,1);
%% collect relevant results
for f = 1 : length(d)
    skip = false;
    fname = d(f).name;
    sfname = split(fname, '_');
    if length(sfname) < 3
        continue
    end
    for v = 1 : n
        if strcmp(sfname{1}, vars{v})
            ffs(v) = ffs(v) + 1;
            break
        end
        if v == n
            skip = true;
        end
    end
    if ~skip
        load([d(f).folder,'/',fname],'mat_devs','mat_mds');
        dev_cell{v, ffs(v)} = mat_devs;
        mds_cell{v, ffs(v)} = mat_mds;
    end
end
% remove empty cells
ff = min(min(ffs, MAX_FILES));
dev_cell(:,ff+1:end) = [];
mds_cell(:,ff+1:end) = [];
% concatenate results into larger matrices
full_devs = cell(n,1);
full_mds  = cell(n,1);
for v = 1 : n
    runs = size(dev_cell{1},3);
    devs = nan(size(dev_cell{1},1),size(dev_cell{1},2),runs*ff);
    mds  = nan(size(mds_cell{1},1),runs*ff);
    for i = 1 : length(dev_cell)
        devs(:,:,(i-1)*runs+1:(i-1)*runs+runs) = dev_cell{v,i};
        mds(:,(i-1)*runs+1:(i-1)*runs+runs) = mds_cell{v,i};
    end
    full_devs{v} = devs;
    full_mds{v}  = mds;
end
%% Compute metrics
% RMSE
rmse_tra = nan(size(full_devs{1},1),n);
rmse_rot = nan(size(full_devs{1},1),n);
rmse_tot = nan(size(full_devs{1},1),n);
for v = 1 : n
    rmse_tra(:,v) = sqrt(mean(mean(full_devs{v}(:,1:3,:).^2,2),3));
    rmse_rot(:,v) = sqrt(mean(mean(full_devs{v}(:,4:6,:).^2,2),3));
    rmse_tot(:,v) = sqrt(mean(mean(full_devs{v}(:,1:6,:).^2,2),3));
    rmse_tra_full = sqrt(mean(rmse_tra(:,v).^2));
    rmse_rot_full = sqrt(mean(rmse_rot(:,v).^2));
    disp(['[I] ',vars{v}, ' RMSE: ', num2str(rmse_tra_full*1e3,'%.2f'),' mm , ',num2str(rad2deg(rmse_rot_full)*1e3,'%.2f'),' mdeg']);
end
% MDS
MCruns = size(dev_cell{1,1},3);
meanMds = nan(size(full_devs{1},1),n);
for v = 1 : n
    meanMds(:,v) = mean(full_mds{v},2);
end
%% Plot results
times = 0 : 1/fp.IMU_Hz : size(full_devs{1},1)/fp.IMU_Hz - 1/fp.IMU_Hz;
% RMSE
figure
% translations
subplot(2,1,1)
plot(times,rmse_tra .* 1e3)
box on
grid on
xlabel('time [s]')
ylabel('RMSE [mm]')
title("Translations")
xlim([times(1), times(end)])
%ylim([0,max(max(rmse_tra) .* 1e3)])
legend(vars,'Location','best')
% orientations
subplot(2,1,2)
plot(times,rad2deg(rmse_rot) .* 1e3)
hold on
box on
grid on
xlabel('time [s]')
ylabel('RMSE [mdeg]')
title("Orientations")
xlim([times(1), times(end)])
%ylim([0,max(max(rad2deg(rmse_rot) .* 1e3))])
legend(vars,'Location','best')
% total pose
figure
plot(times,rmse_tot .* 1e3)
box on
grid on
xlabel('time [s]')
ylabel('RMSE [mm|mdeg]')
title("RMSE Pose (Translations and Orientations)")
xlim([times(1), times(end)])
%ylim([0,max(max(rmse_tra) .* 1e3)])
legend(vars,'Location','best')
%% Mahalanobis distance
figure
plot(times, meanMds)
hold on
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
%ylim([max(min(minMds),1e-0),max([maxMds; mds_limits(2)])])
if max(maxMds) > 10 * mds_limits(2)
    set(gca, 'YScale', 'log')
end
legend(vars,'Location','best')