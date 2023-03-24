function F = buildF_rot(coords, caliLS2IMU)
%BUILDF_ROT jacobian matrix for the var.prop.
% INFO: automatically computed
%   Inputs:
%       - coords: cartesian coordinates in s-frame (nx3)
%       - caliLS2IMU: pose of LiDAR in IMU-frame (6x1,tx,ty,tz,omega,phi,kappa)
%   Output:
%       - F: jacobian matrix (nx3*n)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    x = coords(:,1)';
    y = coords(:,2)';
    z = coords(:,3)';
    
    omLI = caliLS2IMU(4);
    phLI = caliLS2IMU(5);
    kaLI = caliLS2IMU(6);

    tempF = zeros(9, size(coords,1));
    tempF(1,:) = y*(sin(kaLI)*sin(omLI) + cos(kaLI)*cos(omLI)*sin(phLI)) + z*(cos(omLI)*sin(kaLI) - cos(kaLI)*sin(omLI)*sin(phLI));
    tempF(2,:) = z*cos(kaLI)*cos(omLI)*cos(phLI) - x*cos(kaLI)*sin(phLI) + y*cos(kaLI)*cos(phLI)*sin(omLI);
    tempF(3,:) = z*(cos(kaLI)*sin(omLI) - cos(omLI)*sin(kaLI)*sin(phLI)) - y*(cos(kaLI)*cos(omLI) + sin(kaLI)*sin(omLI)*sin(phLI)) - x*cos(phLI)*sin(kaLI);
    tempF(4,:) = -y*(cos(kaLI)*sin(omLI) - cos(omLI)*sin(kaLI)*sin(phLI)) - z*(cos(kaLI)*cos(omLI) + sin(kaLI)*sin(omLI)*sin(phLI));
    tempF(5,:) = z*cos(omLI)*cos(phLI)*sin(kaLI) - x*sin(kaLI)*sin(phLI) + y*cos(phLI)*sin(kaLI)*sin(omLI);
    tempF(6,:) = z*(sin(kaLI)*sin(omLI) + cos(kaLI)*cos(omLI)*sin(phLI)) - y*(cos(omLI)*sin(kaLI) - cos(kaLI)*sin(omLI)*sin(phLI)) + x*cos(kaLI)*cos(phLI);
    tempF(7,:) = y*cos(omLI)*cos(phLI) - z*cos(phLI)*sin(omLI);
    tempF(8,:) = -x*cos(phLI) - z*cos(omLI)*sin(phLI) - y*sin(omLI)*sin(phLI);
    
    colF = reshape(tempF, 3, size(coords, 1) * 3)';

    F = zeros(size(coords, 1) * 3, size(coords, 1) * 3);
    for i = 1 : size(coords, 1)
        F((i-1)*3+1 : (i-1)*3+3, (i-1)*3+1 : (i-1)*3+3) = colF((i-1)*3+1 : (i-1)*3+3,:);
    end
end

