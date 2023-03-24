function [intRegions, pp, fnames, refPts] = importROI(path, tr)
%IMPORTROI imports reference panels from files
%   imports reference panels and sets up variables for there regions, plane
%   parameters, filenames and corner points
%   plane parameters are corrected for possible offsets from the measurement
%   Inputs:
%       - path: path to files (string)
%       - tr: info about the planes in env. (struct)
%   Outputs:
%       - intRegions: volumes of interest around panels (nx6)
%       - pp: plane parameters of reference panels (4xn)
%       - fnames: names of panels (cells -> string)
%       - refPts: corner points (measured?) (cells -> mx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    files = dir(path);
    files = files(~ismember({files.name},{'.','..'}));
    intRegions = nan(length(files), 6); % min_x, max_x, min_y, max_y, min_z, max_z
    if nargout > 1
        pp = nan(4, length(files));
        fnames = cell(length(files),1);
    end
    if nargout > 3
        refPts = cell(length(files),1);
    end
    for i = 1 : length(files)
        % read file
        fname = files(i).name;
        pts = readmatrix([path, fname]);
        pts = pts(:, 2:4) * 1e-3; % mm -> m
        % correct for offset from measurement (CCR or TP tip)
        pts = adjustForOffset(pts, tr.mid_point(1:3), tr.REF_OFFSET);
        refPts{i} = pts;
        intRegions(i,:) = [min(pts(:,1)) - tr.REF_BUFFER, max(pts(:, 1)) + tr.REF_BUFFER, ...
                           min(pts(:,2)) - tr.REF_BUFFER, max(pts(:, 2)) + tr.REF_BUFFER, ...
                           min(pts(:,3)) - tr.REF_BUFFER, max(pts(:, 3)) + tr.REF_BUFFER];
        if nargout > 1
            % compute plane params and standardize normal vector direction
            pp(:, i) = planeEst(pts);
            testDist = pp(:,i)' * tr.mid_point';
            if testDist > 0
                pp(:,i) = pp(:,i) .* -1;
            end
            % save file name for title
            fnames{i} = replace(fname, '_', '\_');
        end
    end
end