% transformations ...
% ... between T-Probe frame t and IMU frame s
H_MS2B = [quaternion(tr.pose_LS_IMU(4:6)','euler','XYZ','point').rotmat('point'),tr.pose_LS_IMU(1:3); 0,0,0,1];
H_MS2B(1:3,4) = H_MS2B(1:3,4) .* 1e-3;
% ... between T-Probe and Body frame
R_TP = rotmat3D(tr.pose_TP(4), tr.pose_TP(5), tr.pose_TP(6));
H_B2TP = [R_TP, tr.pose_TP(1:3)*1e-3; zeros(1,3), 1];

H_MS2TP = H_B2TP * H_MS2B;
ipos = H_MS2TP(1:3,4);
iq = quaternion(H_MS2TP(1:3,1:3),'rotmat','point');

% ... between LT frame l and room frame r
R_WCS = rotmat('z', tr.pose_WCS(6)) * rotmat('y', tr.pose_WCS(5)) * rotmat('x', tr.pose_WCS(4));
H_LT2R  = [R_WCS, tr.pose_WCS(1:3)*1e-3; zeros(1,3), 1];
tr.H_l_r = eye(4);
rot_LT2r = quaternion(H_LT2R(1:3,1:3),'rotmat','point');

% ... between LS frame s and body frame b
R_PCS = rotmat3D(tr.pose_PCS(4), tr.pose_PCS(5), tr.pose_PCS(6));
H_LS2B = [R_PCS, tr.pose_PCS(1:3); zeros(1,3), 1];
% transformation LS pcls from LS to IMU
H_LS2MS = H_MS2B^-1 * H_LS2B;
tr.H_s_t = H_LS2MS;
fp.caliLS2MS = [H_LS2MS(1:3,4)',quaternion(H_LS2MS(1:3,1:3),'rotmat','point').euler('XYZ','point')];

% inital trajectory - in local frame -> transformed in room frame
%  read LT measurements and convert to cartesian coordinates
LT_tracking_obs(:, 2:3) = deg2rad(LT_tracking_obs(:, 2:3));
LT_tracking_obs(:, 2:4) = LT_tracking_obs(:, 4:-1:2);
LT_tracking_obs(:, 2:4) = spherical2cart2(LT_tracking_obs(:, 2:4),true) .* 1e-3; % [mm -> m]

%  transform T-Probe poses into room frame
pos = [LT_tracking_obs(:,2:4)'; ones(1, length(LT_tracking_obs))]; % [mm]
LT_tpos = H_LT2R * pos;
LT_tpos = LT_tpos(1:3,:)'; % [m]
LT_orients = rot_LT2r .* quaternion(LT_tracking_obs(:,5:8));
%  initial values for pose of IMU
H_MS2r = [LT_orients(1,:).rotmat('point'), LT_tpos(1,1:3)'; 0,0,0,1];
ipos_r = H_MS2r * [ipos; 1];
iq_r = LT_orients(1,:) * iq;
% values for the LS
H_LS2r = [iq_r.rotmat('point'), ipos_r(1:3);0,0,0,1] * H_LS2MS;
rot_LS2r = quaternion(H_LS2r(1:3,1:3), 'rotmat', 'point');

initialP = ipos_r(1:3); % [m]
initialO = iq_r; % [-]