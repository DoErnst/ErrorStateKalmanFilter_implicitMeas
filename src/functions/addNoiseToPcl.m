function LS_obs = addNoiseToPcl(LS_obs_simu, st, fp)
%ADDNOISETOPCL adds measurement noise to simulated perfect pcls
%   noise is based on provided information and resolution is applied
%   additionally, a distance offset can be added
%   Inputs:
%       - LS_obs_simu - simulated perfect pcl (t,x,y,z,obj,sc)
%       - st: stochastic model of LiDAR (struct)
%       - fp: functional parameters (struct)
%   Outputs:
%       - LS_obs: LiDAR point cloud with added noise (t,x,y,z)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    LS_obs = nan(size(LS_obs_simu,1),4);
    LS_obs(:,1) = LS_obs_simu(:,1);
    % process in chunks to avoid too large operations for memory
    for i = 1 : 1e6 : size(LS_obs_simu)
        if i + 1e6-1 > size(LS_obs_simu)
            ix_end = size(LS_obs_simu);
        else
            ix_end = i + 1e6-1;
        end
        % measurement elements are spherical coordinates so we transform
        obs = cart2spherical(LS_obs_simu(i:ix_end,2:4));
        % simple noise: normal distributed + distance offset
        noise = [randn(size(obs,1),1) * sqrt(st.VCM_LS(1,1)), ...
                 randn(size(obs,1),1) * sqrt(st.VCM_LS(2,2)), ...
                 randn(size(obs,1),1) * sqrt(st.VCM_LS(3,3))];
        % dist offset: negative of estimated value
        dist_offset = -[(fp.r0_LS + randn()*st.std_r0_LS) * ones(size(obs,1),1), ...
            zeros(size(obs,1),2)];
        obs_n = obs + noise + dist_offset;
        % apply resolution of LS
        %obs_n(:,1) = round(obs_n(:,1).*5,2)./5;
        % transform back to cartesian points for exchange
        LS_obs(i : ix_end, 2:4) = spherical2cart(obs_n);
    end

end

