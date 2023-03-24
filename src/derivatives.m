%% Symbolic for derivations
% stand alone script to determine the partial derivatives of the functional
% model for the filter update and stochastic model
%
% Copyright (c) 2023 Dominik Ernst under MIT License

%% observation level
% transformation spherical -> cartesian
syms tx ty tz di el az vc
t = [tx; ty; tz];
xs = di * cos(el) * cos(az);
ys = di * cos(el) * sin(az);
zs = di * sin(el) + vc;
s2c = [xs; ys; zs];
A_s2c = jacobian(s2c, [di, el, az]);
%% trafo LS -> IMU
syms txLI tyLI tzLI omLI phLI kaLI
RLI = rotmat("z", kaLI) * rotmat("y", phLI) * rotmat("x", omLI);
pLI = [txLI; tyLI; tzLI] + RLI * s2c;
A_vp = jacobian(pLI, [di, el, az]);
% transformed points w/ VCM from propagation before
syms x y z
tf = [txLI; tyLI; tzLI] + RLI * [x;y;z];
A_rot = jacobian(tf, [omLI, phLI, kaLI]);
% transformed points w/ full variance propagation
A_tf = jacobian(tf, [x,y,z, txLI, tyLI, tzLI, omLI, phLI, kaLI]);
%% functional model with rotations using quaternions
syms qw qx qy qz
qv = [qx; qy; qz];
q = [qw; qv];
tfq = t + quat2rotmat(q) * [x; y; z];
% plane params for distance computation
syms ax ay az nd
nvec = [ax; ay; az];
dq = nvec' * tfq + nd;
Afq = jacobian(dq, [tx, ty, tz, qw, qx, qy, qz]); % derived to true state (4x1) q
Bfq = jacobian(dq, [x, y, z]);
A_ext = [Afq(1:3), zeros(1,3), Afq(4:7),zeros(1,6)];
Xdx = blkdiag(eye(6), 0.5 .* [-qv'; qw * eye(3) + axiator(qv)'],eye(6)); % for global error
Aq_full = A_ext * Xdx;

%% functional model with rotations using quaternions (lin. interpolated)
syms tau tx_o ty_o tz_o qw_o qx_o qy_o qz_o
qv_o = [qx_o; qy_o; qz_o];
q_o = [qw_o; qv_o];
% interpolation using tau (assumed to be between 0 (old) and 1 (new))
t_o = [tx_o; ty_o; tz_o]; % new translations
t_k = (1-tau) * t_o + tau * t; % translations for tau
% "classical" SLERP: doesn't work with jacobian
%q_k = quatMult(q_o, rotvec2quat(tau .* quat2rotvec(quatMult(qc(q_o), q))));
% using Sola, 2017: p.29 -> Shoemake (1985)
dth = acos(q_o' * q);
q_k = q_o * (sin((1-tau)*dth) / sin(dth)) + q * (sin(tau*dth)/sin(dth));
tfq_k = t_k + quat2rotmat(q_k) * [x; y; z];
dq_k = nvec' * tfq_k + nd;
Afq_k = jacobian(dq_k, [tx, ty, tz, qw, qx, qy, qz]); % derived to true state (4x1) q
Bfq_k = jacobian(dq_k, [x, y, z]);
A_ext_k = [Afq_k(1:3), zeros(1,3), Afq_k(4:7),zeros(1,6)];
Xdx_k = blkdiag(eye(6), 0.5 .* [-q_k(2:4)'; q_k(1) * eye(3) + axiator(q_k(2:4))'],eye(6)); % for global error
Aq_full_k = A_ext_k * Xdx_k;



%% helper functions

% quaternions
function mat = axiator(vec)
    mat = [0, -vec(3), vec(2);
           vec(3), 0, -vec(1);
           -vec(2), vec(1), 0];
end

function R = quat2rotmat(q)
    qw = q(1);
    qv = q(2:4);
    R = (qw^2 - qv' * qv) * eye(3) + 2 * (qv * qv') + 2 * qw * axiator(qv);
end

% rotation matrices from euler angles
function R = rotmat(axis, angle)
    
    c = cos(angle);
    s = sin(angle);

    if axis == "x" 
        R = [1, 0, 0;
             0, c, -s;
             0, s, c];
    end
    if axis == "y"
        R = [c, 0, s;
             0, 1, 0;
             -s, 0, c];
    end
    if axis == "z"
        R = [c, -s, 0;
             s, c, 0;
             0, 0, 1];
    end
end
function R = rotmat3D(omega, phi, kappa)
    R =  rotmat("z", kappa) * rotmat("y", phi) * rotmat("x", omega);
end

