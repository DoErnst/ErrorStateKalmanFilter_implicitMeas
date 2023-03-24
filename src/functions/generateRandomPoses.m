function random_poses = generateRandomPoses(wp, sizeOfVolume)
%GENERATERANDOMPOSES generates random poses in a given volume
%   generates a given number of random poses in a given room volume
%   generated positions have a minimum distance of 0.5 m to the walls
%   Inputs:
%       - wp: number of positions (scalar)
%       - sizeOfVolume: size of room for generation (x,y,z) [m]
%   Outputs:
%       - random_poses: randomly generated poses (x,y,z,q0,q1,q2,q3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    translations = rand(wp, 3) .* (sizeOfVolume - 2) + 1; % m
    rot_ = rand(wp, 3);
    thetas = acos(2 .* rot_(:,1) - 1);
    phis = 2 .* pi .* rot_(:,2);
    psis = 2 .* pi .* rot_(:,3);
    vecs = [cos(thetas).*cos(phis),cos(thetas).*sin(phis),sin(thetas)];
    qs = [cos(psis./2), sin(psis./2).*vecs];
    qs = qs .* sign(qs(:,1));
    random_poses = [translations, qs];

end