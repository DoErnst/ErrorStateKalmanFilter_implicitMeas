function A = areaTriangle(tri)
%AREATRIANGLE computes the area of a triangle based on the corner points
%   computes the area of a given triangle
%   Inputs:
%       - tri: coordinates of 2D triangle (x,y)
%   Outputs:
%       - A: area of triangle
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    A = abs( ...
            (tri(1,1) * (tri(2,2)-tri(3,2))) + ...
            (tri(2,1) * (tri(3,2)-tri(1,2))) + ...
            (tri(3,1) * (tri(1,2)-tri(2,2))) ...
        ) / 2;
end

