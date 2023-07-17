%% Plots for visualizing results
% automatically parses through output folder and looks for results with the
% corresponding results prefix to visualize the results
%
% Copyright (c) 2023 Dominik Ernst under MIT License

addpath('functions')
importSimulationInputs;
% iterate through result files and collect result matrices
%d = dir([tr.DATAPATH,'/output/']);
d = dir('..\data\output');
dev_cell = cell(1, length(d));
mds_cell = cell(1, length(d));
ff = 0; % counter for found files
for f = 1 : length(d)
    fname = d(f).name;
    sfname = split(fname, '_');
    if length(sfname) < 3 || ~strcmp(sfname{1}, rs.prefix)
        continue
    end
    ff = ff + 1;
    load([d(f).folder,'/',fname],'mat_devs','mat_mds');
    dev_cell{ff} = mat_devs;
    mds_cell{ff} = mat_mds;
end
% remove empty cells
dev_cell(ff+1:end) = [];
mds_cell(ff+1:end) = [];
% concatenate results into larger matrices
runs = size(dev_cell{1},3);
devs = nan(size(dev_cell{1},1),size(dev_cell{1},2),runs*ff);
mds  = nan(size(mds_cell{1},1),runs*ff);
for i = 1 : length(dev_cell)
    devs(:,:,(i-1)*runs+1:(i-1)*runs+runs) = dev_cell{i};
    mds(:,(i-1)*runs+1:(i-1)*runs+runs) = mds_cell{i};
end
%% Plots
times = 0 : 1/fp.IMU_Hz : tr.trajLen;
plotDevsMds(times, devs, mds, rs)