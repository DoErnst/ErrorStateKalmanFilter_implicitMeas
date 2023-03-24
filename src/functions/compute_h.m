function h = compute_h(x_k, coords_k, asgdPlanes, x_k_o, tau)
%BUILDA Computes the misclosure vector
%   INFO: computed by script "derivatives.m"
%   Inputs:
%       - x_k: states of current epoch (16x1)
%       - coords_k: points assigned to panels (nx3) [m]
%       - asgdPlanes: plane parameters of assigned panels (nx4) [-|m]
%       - x_k_o: states of prior epoch (16x1)
%       - tau: vector for interpolation (nx1, between 0 and 1)
%   Outputs:
%       - h: misclosure vector (nx1)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % states (current)
    tx = x_k(1);
    ty = x_k(2);
    tz = x_k(3);
    qw = x_k(7);
    qx = x_k(8);
    qy = x_k(9);
    qz = x_k(10);
    % states (old)
    tx_o = x_k_o(1);
    ty_o = x_k_o(2);
    tz_o = x_k_o(3);
    qw_o = x_k_o(7);
    qx_o = x_k_o(8);
    qy_o = x_k_o(9);
    qz_o = x_k_o(10);
    % observations
    x = coords_k(:,1);
    y = coords_k(:,2);
    z = coords_k(:,3);
    % plane params
    ax = asgdPlanes(:,1);
    ay = asgdPlanes(:,2);
    az = asgdPlanes(:,3);
    nd = asgdPlanes(:,4);

    % check for difference between old and new rotation
    dq = qw.*qw_o + qx.*qx_o + qy.*qy_o + qz.*qz_o;

    if abs(1 - dq) < 1e-10 % avoid numerical problems when no rotation change occurs
        % non-interpolated version
        h = nd + ax.*(tx - y.*(2.*qw.*qz - 2.*qx.*qy) + z.*(2.*qw.*qy + 2.*qx.*qz) + x.*(qw.^2 + qx.^2 - qy.^2 - qz.^2)) + ay.*(ty + x.*(2.*qw.*qz + 2.*qy.*qx) - z.*(2.*qw.*qx - 2.*qy.*qz) - y.*(- qw.^2 + qx.^2 - qy.^2 + qz.^2)) + az.*(tz - x.*(2.*qw.*qy - 2.*qz.*qx) + y.*(2.*qw.*qx + 2.*qz.*qy) - z.*(- qw.^2 + qx.^2 + qy.^2 - qz.^2));
    else
        % non-interpolated version
        h = nd + ax.*(tau.*tx - tx_o.*(tau - 1) + y.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + z.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + x.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2)) + ay.*(tau.*ty - y.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2) - ty_o.*(tau - 1) + x.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + z.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)))) + az.*(tau.*tz - z.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2) - tz_o.*(tau - 1) + x.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + y.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))));
    end

end