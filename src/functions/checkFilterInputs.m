function [st, fp, pm] = checkFilterInputs(st, fp, pm)
%CHECKFILTERINPUTS checks to structs for the filter settings
%   makes sure all settings for the filter are set/plausible
%   throws warning otherwise and sets defaults
%   Inputs:
%       - st: stochastic model for sensors (struct)
%       - fp: functional parameters (trafos, dist offset, etc) (struct)
%       - pm: parameters for the filter (struct)
%       - IMU_obs: accelerations and rotation rates
%           (t,a_x,a_y,a_z,r_x,r_y,r_z)
%   Outputs:
%       - st: corrected stochastic model for sensors (struct)
%       - fp: corrected functional parameters (struct)
%       - pm: corrected parameters for the filter (struct)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % checks for stochastic parameters
    if ~isfield(st, 's_a_w')
        st.s_a_w = [1,1,1];
        warning('s_a_w not set, using default value');
    end
    if ~isfield(st, 's_omega_w')
        st.s_omega_w = [1,1,1];
        warning('s_omega_w not set, using default value');
    end
    if ~isfield(st, 's_a_n')
        st.s_a_n = [1,1,1];
        warning('s_a_n not set, using default value');
    end
    if ~isfield(st, 's_omega_n')
        st.s_omega_n = [1,1,1];
        warning('s_omega_n not set, using default value');
    end
    if ~isfield(st, 'std_r0_LS')
        st.std_r0_LS = 1e-4;
        warning('std_r0_LS not set, using default value');
    end
    if ~isfield(st, 'VCM_LS')
        st.VCM_LS = diag([(0.03^2/3); 1.5e-3^2; 3e-3^2]);
        warning('VCM_LS not set, using default value');
    end
    if ~isfield(st, 'VCM_LS_IMU')
        st.VCM_LS_IMU = zeros(6);
        warning('VCM_LS_IMU not set, using default value');
    end
    % checks for filter settings
    if ~isfield(pm, 'SUBSAMPLING_FACTOR')
        pm.SUBSAMPLING_FACTOR = 1;
    end
    if ~isfield(pm, 'REM_OUTLIERS_LS')
        pm.REM_OUTLIERS_LS = 1;
    end
    % checks for functional parameters
    if ~isfield(fp, 'g_vec')
        fp.g_vec = [0;0;9.81];
        warning('g_vec not set, using default value');
    end
    if ~isfield(fp, 'caliLS2MS')
        fp.caliLS2MS = [0,0,0,0,0,0];
        warning('caliLS2MS not set, using default value');
    end
    if ~isfield(fp, 'r0_LS')
        fp.r0_LS = 0;
    end
    if ~isfield(fp, 'R_corr_IMU')
        fp.R_corr_IMU = eye(3);
    end

end