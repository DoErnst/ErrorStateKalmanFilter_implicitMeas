function refL = plotCoordSys(pos, ori, drawFactor, width, currL)
%PLOTCOORDSYS Plots the origin and orientation of a coordinate system
%   draws a trihedron to show the origin and orientation of a frame in an
%   existing plot or updates them
%   Inputs:
%       - pos: position of the frame (x,y,z)
%       - ori: orientation of the frame (either quaternion or rotation
%           matrix)
%       - drawfactor: length for lines of trihedron (default: 1)
%       - width: line width of lines of trihedron (default: 0.5)
%       - currL: reference to old lines for update (animated plots)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 3
        drawFactor = 1;
    end
    if nargin < 4
        width = 0.5;
    end

    if size(ori,1) < 3
        rm = ori.rotmat('point');
    else
        rm = ori;
    end

    if nargin < 5
        rm = rm .* drawFactor;
        
        refL(1) = line([pos(1),pos(1)+rm(1,1)],[pos(2),pos(2)+rm(2,1)],[pos(3),pos(3)+rm(3,1)],'Color','red','LineWidth',width);
        refL(2) = line([pos(1),pos(1)+rm(1,2)],[pos(2),pos(2)+rm(2,2)],[pos(3),pos(3)+rm(3,2)],'Color','green','LineWidth',width);
        refL(3) = line([pos(1),pos(1)+rm(1,3)],[pos(2),pos(2)+rm(2,3)],[pos(3),pos(3)+rm(3,3)],'Color','blue','LineWidth',width);
    else
        set(currL(1), 'XData', [pos(1),pos(1)+rm(1,1)], 'YData', [pos(2),pos(2)+rm(2,1)], 'ZData', [pos(3),pos(3)+rm(3,1)]);
        set(currL(2), 'XData', [pos(1),pos(1)+rm(1,2)], 'YData', [pos(2),pos(2)+rm(2,2)], 'ZData', [pos(3),pos(3)+rm(3,2)]);
        set(currL(3), 'XData', [pos(1),pos(1)+rm(1,3)], 'YData', [pos(2),pos(2)+rm(2,3)], 'ZData', [pos(3),pos(3)+rm(3,3)]);
    end

end

