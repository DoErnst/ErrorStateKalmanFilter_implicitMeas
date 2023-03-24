function [intRegions, pp, H_r_p, refPts_r] = importROI2(refPts, tr)
%IMPORTROI2 imports reference panels 
%   imports reference panels and sets up variables for there regions, plane
%   parameters, filenames and corner points
%   plane parameters are corrected for possible offsets from the measurement
%   Inputs:
%       - refPts: corner points of panels (cells -> mx3)
%       - tr: info about the planes in env. (struct)
%   Outputs:
%       - intRegions: volumes of interest around panels (nx6)
%       - pp: plane parameters of reference panels (4xn)
%       - fnames: names of panels (cells -> string)
%       - refPts_r: corner points in plane-frame (cells -> mx4)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    intRegions = nan(length(refPts), 6); % min_x, max_x, min_y, max_y, min_z, max_z
    if nargout > 1
        pp = nan(4, length(refPts));
    end
    for i = 1 : length(refPts)
        % read file
        pts = refPts{i};
        % correct for offset from measurement (CCR or TP tip)
        pts = adjustForOffset(pts, tr.mid_point(1:3), tr.REF_OFFSET);
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
        end
    end
    n_r = size(tr.refPts,1);
    H_r_p = zeros(4,4,n_r);
    refPts_r = cell(n_r,1);
    for i = 1 : n_r
        H_r_p(1:3,4,i) = mean(tr.refPts{i});
        [~,~,H_r_p(1:3,1:3,i)] = planeEst(tr.refPts{i});
        H_r_p(4,4,i) = 1;
        H_r_p(:,:,i) = H_r_p(:,:,i)^-1;
        refPts_r{i,1} = [tr.refPts{i}, ones(size(tr.refPts{i},1),1)] * H_r_p(:,:,i)';
    end
end