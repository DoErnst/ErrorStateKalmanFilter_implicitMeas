function [coords_ix, asgdPlanes, p_ix] = findPtsInROIforFilter(tr, rpcl, remOutliers)
%FINDPTSINROIFORFILTER finds points in ROI for assignment
%   steps for assignment of points to planes in env:
%       - computation of distances between all points and planes
%       - set up matrix showing which points are in volumes of which ref
%           panels
%       - filter out distances if point isn't in volume of ref panel
%       - find minimum distances from points to planes
%       - filter out points without assignments
%       - per ref panel:
%           - find all points assigned to it
%           - perform outlier search using median absolute deviation
%   Info: for efficient computation in larger environments, the
%       environmental info should be preselected/ reduced using spatial
%       index structures
%   Inputs:
%       - tr: info about the planes in env. (struct)
%       - rpcl: point cloud in r-frame (same as planes) (nx3, (x,y,z))
%       - intPosLiDAR: origin of rays (sensor position) (nx3, (x,y,z))
%       - remOutliers: should outliers be removed (logical)
%   Outputs:
%       - coords_ix: indices of points assigned (nx1)
%       - asgdPlanes: planes parameters for assigned points (nx4)
%       - p_ix: indices of assigned to planes (nx1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    pcl = rpcl(1:3,:)';
    % collect indices of points and the corresponding pp in cell arrays
    tempC = cell(size(tr.intRegions, 1),1);
    tempP = cell(size(tr.intRegions, 1),1);
    tempI = cell(size(tr.intRegions, 1),1);
    % compute distances -> NECESSARY??
    distsMat = rpcl' * tr.intPP;
    % find correspondences: in regions
    logIntRegionMat = false(size(pcl,1), size(tr.intPP,2));
    for r = 1 : size(tr.intRegions, 1)
        logIntRegionMat(:,r) = pointsInRegion(pcl, tr.intRegions(r, :));
    end
    % ignore all points outside the ROIs
    distsMat(~logIntRegionMat) = NaN;
    % get indices of minimum distances per point
    [~,ix] = min(abs(distsMat),[],2);
    % filter point without valid assignment
    d_ = sum(isnan(distsMat),2);
    ix(d_ == size(tr.intPP,2)) = NaN;
    %   combine both steps and only refpanel with minimum distance
    % afterwards outlier search as implemented
    f_ix = unique(ix);
    for r = f_ix'
        if isnan(r)
            continue
        end
        temp_ix = find(ix == r); % get indices of pts in ROI
        % check for inliers by computing distances of pts in ROI to pp
        if remOutliers
            regionPcl = [pcl(ix == r,:), ones(sum(ix == r),1)]';
            dists = tr.intPP(:,r)' * regionPcl;
            logDist = ~isoutlier(dists, 'median')'; % mark inliers
            tempC{r} = temp_ix(logDist)'; % only use ix of inliers
            tempP{r} = repmat(tr.intPP(:,r),1,sum(logDist));
            tempI{r} = repmat(r,1,sum(logDist));
        else
            tempC{r} = temp_ix';
            tempP{r} = repmat(tr.intPP(:,r),1,sum(logIntRegion));
            tempI{r} = repmat(r,1,sum(logIntRegion));
        end
    end
    % concatenate all point with their corresponding planes
    coords_ix = [tempC{:}]';
    asgdPlanes = [tempP{:}]';
    p_ix = [tempI{:}]';
end