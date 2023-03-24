function traj = generateSimpleTraj(type)
%GENERATESIMPLETRAJ generates a simple 6-DoF trajectory
%   uses waypoints to generate simple trajectories (currently 3 types)
%   Inputs:
%       - type: switches between different (hard-coded) trajs (scalar)
%   Outputs:
%       - traj: trajectory (t,x,y,z,q0,q1,q2,q3)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    switch type
        case 1 % movement with simultaneous rotations
            transl = [1,1,1.5; % Wp 1
                      1,1,1.5; % WP 1 repeated (for standing start)
                      9,1,1.5; % WP 2
                      9,7,1.5; % WP 3
                      1,7,1.5; % WP 4
                      1,1,1.5; % WP 1 again (full round)
                      1,1,1.5]; % WP 1 again again (standing stop)
            rots   = [1,0,0,0;
                      1,0,0,0;
                      sqrt(2)/2,0,0,sqrt(2)/2;
                      0,0,0,1;
                      -sqrt(2)/2,0,0,sqrt(2)/2;
                      -1,0,0,0;
                      -1,0,0,0];
            times  = [0; 5; 13; 19; 27; 33; 38];
        case 2 % movement and rotations alternatingly
            transl = [1,1,1.5; % WP 1
                      1,1,1.5; % WP 1 repeated (for standing start)
                      9,1,1.5; % WP 2
                      9,1,1.5; % WP 2 (for rotation)
                      9,7,1.5; % WP 3
                      9,7,1.5; % WP 3 (for rotation)
                      1,7,1.5; % WP 4
                      1,7,1.5; % WP 4 (for rotation)
                      1,1,1.5; % WP 1 again (full round)
                      1,1,1.5; % WP 1 again (for rotation)
                      1,1,1.5; % WP 1 again (for rotation)
                      1,1,1.5; % WP 1 again again (standing stop)
                      1,1,2.5]; % vertical movement
            rots   = [1,0,0,0;
                      1,0,0,0;
                      1,0,0,0;
                      sqrt(2)/2,0,0,sqrt(2)/2;
                      sqrt(2)/2,0,0,sqrt(2)/2;
                      0,0,0,1;
                      0,0,0,1;
                      -sqrt(2)/2,0,0,sqrt(2)/2;
                      -sqrt(2)/2,0,0,sqrt(2)/2;
                      -1,0,0,0;
                      -sqrt(2),0,sqrt(2),0;
                      -sqrt(2),0,sqrt(2),0;
                      -sqrt(2),0,sqrt(2),0];
            times  = [0; 5; 13; 18; 24; 29; 37; 42; 48; 53; 58; 60; 65];
        otherwise % stationary
            transl = [1,1,1.5; % Wp 1
                      1,1,1.5];
            rots   = [1,0,0,0;
                      1,0,0,0];
            times  = [0; 10];
    end
    traj = [times, transl, rots];


end

