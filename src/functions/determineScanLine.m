function abs_line = determineScanLine(omegas)
%DETERMINESCANLINE determines the scan lines for a velodyne scan
%   for configuration of Velodyne Puck (VLP-16)
%   Input:
%	- omegas: elevation angles of points (nx1)
%   Output:
%	- abs_line: absolute ring number per point (nx1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    angles = deg2rad([1,3,5,7,9,11,13,15]);
    diffs = abs(omegas) - angles; % column-row -> mat
    [~, abs_line] = min(abs(diffs), [], 2); % -> indices of min values
end

