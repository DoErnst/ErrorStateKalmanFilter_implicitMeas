function inTri = inTriangleAng(p, tri)
% INTRIANGLEANG checks if a point lies in a triangle
%   computes the angles between sides of original triangle and points
%   sum of angles should be 2pi otherwise point outside
%   Note: somehow slower than area comparison, still included
%   Inputs:
%       - p: point (2x1)
%       - tri: triangle described by 3 corner points (3x2)
%   Outputs:
%       - inTri: true if p is in triangle (boolean)
%
% Copyright (c) 2023 Dominik Ernst under MIT License


    n1 = (tri(1,:) - p) / norm(tri(1,:) - p);
    n2 = (tri(2,:) - p) / norm(tri(2,:) - p);
    n3 = (tri(3,:) - p) / norm(tri(3,:) - p);

    ang1 = atan2(norm(det([n2;n1])), dot(n1,n2));
    ang2 = atan2(norm(det([n3;n2])), dot(n2,n3));
    ang3 = atan2(norm(det([n1;n3])), dot(n3,n1));

    inTri = abs(mod(ang1 + ang2 + ang3, 2*pi)) < 1e-10;

end