function refPts = generateEnv(dimX,dimY,dimZ)
%GENERATEENV generates an environment of panels
%   generates an environment with panels distributed for further simulations
%   this function is semi-flexible wrt. the passed parameters, use with care
%   Inputs:
%       - dimX: dimension along X (scalar) [m]
%       - dimY: dimension along Y (scalar) [m]
%       - dimZ: dimension along Z (scalar) [m]
%   Outputs:
%       - refPts: matrix of points on reference panels (cell:66x1 -> 4x3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 1
        dimX = 10;
        dimY = 8;
        dimZ = 6;
    end
    
    refPts = cell(66,1);
    % boundaries/ outer walls
    refPts{1} = [0,0,0;dimX,0,0;dimX,dimY,0;0,dimY,0];
    refPts{2} = [0,0,0;0,dimY,0;0,dimY,dimZ;0,0,dimZ];
    refPts{3} = [0,0,0;0,0,dimZ;dimX,0,dimZ;dimX,0,0];
    refPts{4} = [0,0,dimZ;dimX,0,dimZ;dimX,dimY,dimZ;0,dimY,dimZ];
    refPts{5} = [dimX,0,0;dimX,dimY,0;dimX,dimY,dimZ;dimX,0,dimZ];
    refPts{6} = [0,dimY,0;0,dimY,dimZ;dimX,dimY,dimZ;dimX,dimY,0];
    % wall structure: X/Y
    refPts{7} = [1,3,0;2,3,.5;2,5,.5;1,5,0];
    refPts{8} = [2,2,0;4,2,0;4,3,.5;2,3,.5];
    refPts{9} = mirror(refPts{7},[3,0,0]);
    refPts{10} = mirror(refPts{8},[0,4,0]);
    refPts{11} = [2,3,.5;4,3,.5;4,5,.5;2,5,.5];
    refPts(12:16) = mirror_mult(refPts(7:11), [dimX/2,0,0]);
    refPts(17:26) = mirror_mult(refPts(7:16), [0,0,dimZ/2]);
    % wall structure: Y/Z
    refPts{27} = [0,1,2.5;.5,2,2.5;.5,2,3.5;0,1,3.5];
    refPts{28} = [0,2,1.5;0,3,1.5;.5,3,2.5;.5,2,2.5];
    refPts{29} = mirror(refPts{27},[0,2.5,0]);
    refPts{30} = mirror(refPts{28},[0,0,3]);
    refPts{31} = [.5,2,2.5;.5,3,2.5;.5,3,3.5;.5,2,3.5];
    refPts(32:36) = mirror_mult(refPts(27:31), [0,dimY/2,0]);
    refPts(37:46) = mirror_mult(refPts(27:36), [dimX/2,0,0]);
    % wall structure: X/Z
    refPts{47} = [1,0,2.5;2,.5,2.5;2,.5,3.5;1,0,3.5];
    refPts{48} = [2,0,1.5;4,0,1.5;4,.5,2.5;2,.5,2.5];
    refPts{49} = mirror(refPts{47},[3,0,0]);
    refPts{50} = mirror(refPts{48},[0,0,3]);
    refPts{51} = [2,.5,2.5;4,.5,2.5;4,.5,3.5;2,.5,3.5];
    refPts(52:56) = mirror_mult(refPts(47:51), [dimX/2,0,0]);
    refPts(57:66) = mirror_mult(refPts(47:56), [0,dimY/2,0]);

end
%% HELPER FUNCTIONS

function c_pts = mirror_mult(pts_c, axis)
    c_pts = cell(1,length(pts_c));
    for i = 1 : length(pts_c)
        c_pts{i} = mirror(pts_c{i}, axis);
    end
end

function rpts = mirror(pts, axis)
    spts = pts - axis;
    spts(:,axis > 0) = -spts(:,axis > 0);
    rpts = spts + axis;
end