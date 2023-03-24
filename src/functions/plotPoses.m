function plotPoses(poses, logT)
%PLOTPOSES plots the poses of a trajectory in subplots
%   Note: orientations are euler angles
%   Inputs:
%       - poses: trajectory with poses over time (nx6)
%       - logT: marks specific epochs (nx1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin > 1
        n = size(poses, 2);
        colors = zeros(n, 3);
        colors(logT, 2) = 1;
        colors(~logT, 1) = 1;
    end
    figure
    tiledlayout(2,3, 'TileSpacing','compact','Padding', 'compact');
    y_str = {'x [mm]','y [mm]','z [mm]', '\omega [°]', '\phi [°]', '\kappa [°]'};
    for i = 1 : 6
        nexttile
        if i < 4
            plot((1:length(poses))/50,poses(i,:))
        else
            plot((1:length(poses))/50,rad2deg(poses(i,:)))
        end
        xlabel('time since start [s]')
        ylabel(y_str{i})
        xlim([1, length(poses)/50])
        grid on
        if nargin > 1
            hold on
            yl = ylim;
            scatter([1:n]./50, yl(2), 10, colors, '.');
            ylim(yl)
        end
    end

end

