function pts = adjustForOffset(pts, pos, offset)
%ADJUSTFOROFFSET adjusts for the offset in point measurement for planes
%   Inputs:
%       - pts: list of points belonging to one plane (nx3)
%       - pos: central position from which the offset should be applied (3x1)
%       - offset: value of the offset (scalar)
%   Outputs:
%       - pts: list of corrected (shifted) points (nx3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % make sure pos is a column vector
    pos = reshape(pos, 3, 1);
    
    % estimate plane parameters
    pp = planeEst(pts);
    % compute reference distance and test distance
    distRef = norm(pts(1,:)' - pos);
    distTest = norm((pts(1,:)' + pp(1:3) * offset) - pos);
    % adjust coordinates along normal vector of plane
    pts = pts + pp(1:3)' * offset * sign(distTest - distRef);
end

