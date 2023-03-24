function inPoly = inPoly(p, poly)
% INPOLY checks if a point lies in a polygon
%   currently only works for simple (convex) polygons with 4 points
%   strategy: split the 4 points into 2 triangles
%           check if the point is in either triangle
%   NOTE: care for order of points describing polygon (loops!)
%   Inputs:
%       - p: point (2x1)
%       - poly: polygon described by 4 corner points (4x2)
%   Outputs:
%       - inPoly: true if p is in poly (boolean)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % both inputs are assumed to be row vectors
    if size(p,2) < 2
        p = p';
    end
    if size(poly, 2) > 2
        poly = poly';
    end

    inPoly = inTriangle(p, poly(1:3,:)) || inTriangle(p, poly([1,4,3], :));
    
end

