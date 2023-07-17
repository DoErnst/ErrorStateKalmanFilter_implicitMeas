%% Simulations -> visualized
% script for single runs of the filter to visualize and analyze results
%
% Copyright (c) 2023 Dominik Ernst under MIT License

addpath('functions')
importSimulationInputs; % load settings
plottol = 0.5;
FPS = 30;
PLAYBACK_SPEED = 1; % 1 = normal
saveVideo = false;
%% generate trajectory and sensor data
rng(5)
traj = generateRandomTraj(tr.trajLen, [10,8,6]);
%traj = generateSimpleTraj(2);
[IMU_obs, LS_obs, traj_IMU] = generateSensorDataFromTraj(traj, tr, fp);
% add noise to observations
LS_obs_n = addNoiseToPcl(LS_obs, st, fp);
IMU_obs_n = addNoiseToIMU(IMU_obs, traj_IMU, st, fp);
%% filter
% initialize filter
%   compute accelerometer bias (rotate to eliminate grav influence)
R_acc = quaternion(traj_IMU(1,5:8)).rotmat('point');
acc_b = mean(R_acc * IMU_obs_n(1:300, 2:4)',2) - fp.g_vec;
%   setup initial state vector
x0 = [traj(1,2:4), zeros(1,3), traj(1,5:8), acc_b', mean(IMU_obs_n(1:300, 5:7))]';
Qx0 = diag([0.0001^2 * ones(1,3), 0.001^2 * ones(1,3), deg2rad(0.01)^2 * ones(1,3), ...
       st.s_a_w.^2, st.s_omega_w.^2]);
% run filter
LS_obs_n = [LS_obs_n, LS_obs(:,5)];
[states, statesVCM, time] = ESKF_iM_LiDAR_IMU(LS_obs_n,IMU_obs_n,x0, Qx0, st, fp, pm, tr);
% compute deviations between true and filter trajectory
[devs, mds, mds_limits] = compareTraj(traj_IMU(:,2:end),states, statesVCM);
%% Static Plots
plotDevsMds(traj_IMU(:,1),devs,mds);
return % remove for animated plot
%% Animated Plot
% timesteps of LiDAR columns (max. 16 points at once)
timesteps = unique(LS_obs(:,1));
% interpolate true (IMU) traj
traj_int = interp1(traj_IMU(:,1),traj_IMU(:,2:8), timesteps,'pchip');
traj_int(:, 4:7) = traj_int(:, 4:7) ./ vecnorm(traj_int(:,4:7),2,2);
traj_int = [timesteps, traj_int];
% interpolate filter traj
traj_int_f = interp1(traj_IMU(:,1),[states(:,1:3),states(:,7:10)], timesteps,'pchip');
traj_int_f(:, 4:7) = traj_int_f(:, 4:7) ./ vecnorm(traj_int_f(:,4:7),2,2);
traj_int_f = [timesteps, traj_int_f];
LS_obs_trafo_f = nan(size(LS_obs,1), 4);
% for plotting
curr_ix = 1;                                            
old_ix = 1;
ix_old_plot = 1;
len_buffer = ceil(FPS/10) + 1;
circ_buffer_pcl = cell(1, len_buffer);
curr_pcl_ix = 1;
figure
tic;
t_last = toc;
if saveVideo
    vw = VideoWriter('ani','MPEG-4');
    vw.FrameRate = FPS;
    open(vw);
end
for i = 1 : size(traj_int, 1)
    H_IMU = [quat2rotmat(traj_int(i,5:8)),traj_int(i,2:4)';0,0,0,1];
    H_f = [quat2rotmat(traj_int_f(i,5:8)),traj_int_f(i,2:4)';0,0,0,1];
    while LS_obs(curr_ix,1) < timesteps(i)
        curr_ix = curr_ix + 1;
    end
    pts = [LS_obs(old_ix:curr_ix,2:4), ones(curr_ix - old_ix+1,1)];
    LS_obs_trafo_f(old_ix:curr_ix,:) = pts * H_f';
    % for plotting
    if LS_obs(curr_ix,1) > LS_obs(ix_old_plot) + PLAYBACK_SPEED / FPS
        t_now = toc;
        pause(PLAYBACK_SPEED / FPS - (t_now - t_last));
        pix = ix_old_plot : curr_ix;
        circ_buffer_pcl{curr_pcl_ix} = LS_obs_trafo_f(pix, :)';
        c_p_ix = find(traj_IMU(:,1) > LS_obs(curr_ix,1),1,"first");
        if old_ix == 1
            % plot point cloud(s)
            subplot(3,2,[1,5])
            ppcl = plot3(circ_buffer_pcl{1}(1,:),circ_buffer_pcl{1}(2,:),circ_buffer_pcl{1}(3,:),'b.');
            hold on
            filtCS = plotCoordSys(traj_int_f(i,2:4), quaternion(traj_int_f(i,5:8)),0.5,2);
            % true pose
            trueCS = plotCoordSys(traj_int(i,2:4), quaternion(traj_int(i,5:8)));
            plotPanels(tr)
            hold off
            axis equal
            grid on
            box on
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            xlim([0-plottol,10+plottol])
            ylim([0-plottol,8+plottol])
            zlim([0-plottol,6+plottol])
            titlePcl = title(['Filter Trajectory: Time ', num2str(LS_obs(curr_ix,1), '%.2f'),' s']);
            % deviation plots
            % translation
            subplot(3,2,2)
            % complete lines
            plot(traj_IMU(:,1),devs(:,1),'Color',[1,0.8,0.8],'LineWidth',0.1)
            hold on
            plot(traj_IMU(:,1),devs(:,2),'Color',[0.8,1,0.8],'LineWidth',0.1)
            plot(traj_IMU(:,1),devs(:,3),'Color',[0.8,0.8,1],'LineWidth',0.1)
            % current lines
            pcurr(1) = plot(traj_IMU(1:c_p_ix,1),devs(1:c_p_ix,1),'r','LineWidth',1);
            pcurr(2) = plot(traj_IMU(1:c_p_ix,1),devs(1:c_p_ix,2),'g','LineWidth',1);
            pcurr(3) = plot(traj_IMU(1:c_p_ix,1),devs(1:c_p_ix,3),'b','LineWidth',1);
            % progress line
            pline(1) = line([traj_IMU(c_p_ix,1),traj_IMU(c_p_ix,1)],[-10,10],'Color',[0.8,0.8,0.8]);
            hold off
            box on
            grid on
            xlabel('time [s]')
            ylabel('deviation [m]')
            title("Translations")
            xlim([traj_IMU(1,1), traj_IMU(end,1)])
            ylim([max(min(min(devs(:,1:3))),-2),min(max(max(devs(:,1:3))),2)])
            % orientation
            subplot(3,2,4)
            % complete lines
            plot(traj_IMU(:,1),rad2deg(devs(:,4)),'Color',[1,0.8,0.8],'LineWidth',0.1)
            hold on
            plot(traj_IMU(:,1),rad2deg(devs(:,5)),'Color',[0.8,1,0.8],'LineWidth',0.1)
            plot(traj_IMU(:,1),rad2deg(devs(:,6)),'Color',[0.8,0.8,1],'LineWidth',0.1)
            % current lines
            pcurr(4) = plot(traj_IMU(1:c_p_ix,1),rad2deg(devs(1:c_p_ix,4)),'r','LineWidth',1);
            pcurr(5) = plot(traj_IMU(1:c_p_ix,1),rad2deg(devs(1:c_p_ix,5)),'g','LineWidth',1);
            pcurr(6) = plot(traj_IMU(1:c_p_ix,1),rad2deg(devs(1:c_p_ix,6)),'b','LineWidth',1);
            % progress line
            pline(2) = line([traj_IMU(c_p_ix,1),traj_IMU(c_p_ix,1)],[-10,10],'Color',[0.8,0.8,0.8]);
            hold off
            box on
            grid on
            xlabel('time [s]')
            ylabel('deviation [°]')
            title("Orientations")
            xlim([traj_IMU(1,1), traj_IMU(end,1)])
            ylim([max(rad2deg(min(min(devs(:,4:6)))),-10),min(rad2deg(max(max(devs(:,4:6)))),10)])
            % Mahalanobis distance
            subplot(3,2,6)
            % complete lines
            plot(traj_IMU(:,1),mds(:,1),'Color',[0.8,0.8,1],'LineWidth',0.1)
            hold on
            % current lines
            pcurr(7) = plot(traj_IMU(1:c_p_ix,1),mds(1:c_p_ix,1),'b','LineWidth',1);
            % progress line
            pline(3) = line([traj_IMU(c_p_ix,1),traj_IMU(c_p_ix,1)],[1e-10,max(mds)],'Color',[0.8,0.8,0.8]);
            % interval lines
            line([traj_IMU(1,1), traj_IMU(end,1)],[mds_limits(1),mds_limits(1)],'Color',[0,1,0])
            line([traj_IMU(1,1), traj_IMU(end,1)],[mds_limits(2),mds_limits(2)],'Color',[0,1,0])
            hold off
            box on
            grid on
            xlabel('time [s]')
            ylabel('sqrt(d_m) [-]')
            title("sqrt of Mahalanobis distance")
            xlim([traj_IMU(1,1), traj_IMU(end,1)])
            ylim([max(min(mds),1e-0),max(mds)])
            if max(mds) > 10 * mds_limits(2)
                set(gca, 'YScale', 'log')
            end
            if saveVideo
                set(gcf, 'units','pixel','outerposition',[0 0 1600 1200]);
            end
        else % Update data of plots
            % 3D plot of point clouds
            currPcl = [circ_buffer_pcl{:}];
            set(ppcl, 'XData', currPcl(1,:), 'YData', currPcl(2,:), 'ZData', currPcl(3,:));
            plotCoordSys(traj_int_f(i,2:4), quaternion(traj_int_f(i,5:8)),0.5,2,filtCS);
            plotCoordSys(traj_int(i,2:4), quaternion(traj_int(i,5:8)),1,0.5,trueCS);
            titlePcl.String = ['Filter Trajectory: Time ', num2str(LS_obs(curr_ix,1), '%.2f'),' s'];
            % line plots with devs and mds
            set(pcurr(1), 'XData', traj_IMU(1:c_p_ix,1), 'YData', devs(1:c_p_ix,1));
            set(pcurr(2), 'XData', traj_IMU(1:c_p_ix,1), 'YData', devs(1:c_p_ix,2));
            set(pcurr(3), 'XData', traj_IMU(1:c_p_ix,1), 'YData', devs(1:c_p_ix,3));
            set(pcurr(4), 'XData', traj_IMU(1:c_p_ix,1), 'YData', rad2deg(devs(1:c_p_ix,4)));
            set(pcurr(5), 'XData', traj_IMU(1:c_p_ix,1), 'YData', rad2deg(devs(1:c_p_ix,5)));
            set(pcurr(6), 'XData', traj_IMU(1:c_p_ix,1), 'YData', rad2deg(devs(1:c_p_ix,6)));
            set(pcurr(7), 'XData', traj_IMU(1:c_p_ix,1), 'YData', mds(1:c_p_ix,1));
            set(pline(1), 'XData', [traj_IMU(c_p_ix,1), traj_IMU(c_p_ix,1)], 'YData', [-10,10]);
            set(pline(2), 'XData', [traj_IMU(c_p_ix,1), traj_IMU(c_p_ix,1)], 'YData', [-10,10]);
            set(pline(3), 'XData', [traj_IMU(c_p_ix,1), traj_IMU(c_p_ix,1)], 'YData', [1e-10,10000]);
        end
        frame = getframe(gcf);
        if saveVideo
            writeVideo(vw, frame);
        end
        ix_old_plot = curr_ix;
        old_ix = curr_ix + 1;
        curr_pcl_ix = mod(curr_pcl_ix, len_buffer) + 1;
        t_last = t_now;
    end
end
if saveVideo
    close(vw)
end