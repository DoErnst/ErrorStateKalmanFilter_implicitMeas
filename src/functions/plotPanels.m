function plotPanels(tr)
%PLOTPANELS plots panel as patches
%   Inputs:
%       - tr: info about the planes in env. (struct)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    for i = 1 : length(tr.refPts)
        patch(tr.refPts{i}(:,1),tr.refPts{i}(:,2),tr.refPts{i}(:,3),'k', ...
            'FaceAlpha', 0.2)
    end

end