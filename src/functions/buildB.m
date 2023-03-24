function B = buildB(x_k, coords_k, asgdPlanes, x_k_o, tau)
%BUILDB Sets up the restriction matrix
%   INFO: derivatives computed by script "derivatives.m"
%   Inputs:
%       - x_k: states of current epoch (16x1)
%       - coords_k: points assigned to panels (nx3) [m]
%       - asgdPlanes: plane parameters of assigned panels (nx4) [-|m]
%       - x_k_o: states of prior epoch (16x1)
%       - tau: vector for interpolation (nx1, between 0 and 1)
%   Outputs:
%       - B: restriction matrix (nx3n)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % states (current)
    qw = x_k(7);
    qx = x_k(8);
    qy = x_k(9);
    qz = x_k(10);
    % states (old)
    qw_o = x_k_o(7);
    qx_o = x_k_o(8);
    qy_o = x_k_o(9);
    qz_o = x_k_o(10);
    % plane params
    ax = asgdPlanes(:,1);
    ay = asgdPlanes(:,2);
    az = asgdPlanes(:,3);

    tempB = nan(size(coords_k));

    % check for difference between old and new rotation
    dq = qw.*qw_o + qx.*qx_o + qy.*qy_o + qz.*qz_o;

    if abs(1 - dq) < 1e-10 % avoid numerical problems when no rotation change occurs
        % non-interpolated version
        tempB(:,1) = ax.*(qw.^2 + qx.*qx - qy.*qy - qz.*qz) + ay.*(2.*qw.*qz + 2.*qy.*qx) - az.*(2.*qw.*qy - 2.*qz.*qx);
        tempB(:,2) = az.*(2.*qw.*qx + 2.*qz.*qy) - ax.*(2.*qw.*qz - 2.*qx.*qy) - ay.*(- qw.^2 + qx.*qx - qy.*qy + qz.*qz);
        tempB(:,3) = ax.*(2.*qw.*qy + 2.*qx.*qz) - ay.*(2.*qw.*qx - 2.*qy.*qz) - az.*(- qw.^2 + qx.*qx + qy.*qy - qz.*qz);
    else
        % interpolated version
        tempB(:,1) = ay.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + az.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + ax.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2);
        tempB(:,2) = ax.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + az.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) - ay.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2);
        tempB(:,3) = ax.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) + ay.*(2.*((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((2.*qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (2.*qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2))) - az.*(((sin(conj(acos(dq)).*(tau - 1)).*qx_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qx)/conj((1 - (dq).^2).^(1/2))).*((qx_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qx.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) + ((sin(conj(acos(dq)).*(tau - 1)).*qy_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qy)/conj((1 - (dq).^2).^(1/2))).*((qy_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qy.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((sin(conj(acos(dq)).*(tau - 1)).*qz_o)/conj((1 - (dq).^2).^(1/2)) - (sin(conj(acos(dq)).*tau).*qz)/conj((1 - (dq).^2).^(1/2))).*((qz_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qz.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)) - ((qw_o.*sin(acos(dq).*(tau - 1)))/(1 - (dq).^2).^(1/2) - (qw.*sin(tau.*acos(dq)))/(1 - (dq).^2).^(1/2)).^2);
    end

    %% setup up block diagonal matrix
    B = zeros(size(coords_k,1), size(coords_k,1)*3);
    for i = 1 : size(B)
        B(i, (i-1)*3+1 : (i-1)*3+3) = tempB(i,:);
    end
    
end