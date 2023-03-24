function inTri = inTriangle(p, tri, A)
% INTRIANGLE checks if a point lies in a triangle
%   computes the area for the triangles between the original sides and the
%   point and compares to the whole area (can be passed as A)
%   Inputs:
%       - p: point (2x1)
%       - tri: triangle described by 3 corner points (3x2)
%       - A: area of original triangle (scalar)
%   Outputs:
%       - inTri: true if p is in triangle (boolean)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 3
        A  = areaTriangle(tri);
    end
    A1 = areaTriangle([p; tri(2:3,:)]);
    A2 = areaTriangle([p; tri(1:2,:)]);
    A3 = areaTriangle([p; tri([1,3],:)]);

    inTri = A >= A1 + A2 + A3 - 1e-10; % w/ epsilon for floating point err.

end

