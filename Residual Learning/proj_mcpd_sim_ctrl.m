%% Control script for proj_mcpd_sim
clear
close all
reset(gpuDevice)

%% Load Necessary Files
% load("net_AWGN")
load("net_AWGN_Rayleigh")
dnn_network = net.net;

%% Inputs
dist_mat = [10*ones(160,1),linspace(10.05,18,160).'];
att_mat = repmat([1,sqrt(10^(-3/10))],size(dist_mat,1),1);

SNR = 15;                   % SNR to simulate in dB
num_runs = 20;              % Sets the number of trials to do for each distance
isRayleigh = 1;             % Specify whether or not to ignore the att_mat and generate a Rayleigh fading channel.

%% Processing one distance at a time
clear results_all stats_all cal_dist_theo

tic
ctr = 1;
for k = 1:size(dist_mat,1)
    % Run the ranging script
    [results_all(ctr),stats_all(ctr)] =...
        proj_mcpd_sim(num_runs,dist_mat(k,:),att_mat(k,:),...
        SNR,dnn_network,S,isRayleigh);
    ctr = ctr + 1;
    clc
    disp(ctr-1 + "/" + size(dist_mat,1));
end
toc

%% Plot Statistics
plot_statistics(results_all,stats_all,dist_mat,SNR);

%% Functions
function plot_statistics(results_all,stats_all,dist_mat,SNR)
%% Concatenate structures
results = cont_struct(results_all);
stats = cont_struct(stats_all);

% Define a plot name

plot_filename = strcat("Simulated: SNR = ",num2str(SNR)," dB");

% Plot Stats
figure(1);
clf(1)
% Estimated distances
subplot(2,2,1);
hold on
scatter(dist_mat(:,2),100*stats.MUSIC_mean,'o')
scatter(dist_mat(:,2),100*stats.MUSIC_DNN_mean,'x')
scatter(dist_mat(:,2),100*(dist_mat(1,1)),'k.')
hold off
grid on
grid minor
title('Distance Estimations','interpreter','latex')
xlabel('Multi-Path Distance (m)','interpreter','latex')
ylabel('Estimated Distance (cm)','interpreter','latex')
legend({'MUSIC',...
    'MUSIC-DNN',...
    'Theoretical'},...
    'interpreter','latex',...
    'location','best');

% Standard Deviations
subplot(2,2,2);
hold on
scatter(dist_mat(:,2),100*stats.MUSIC_STD,'o')
scatter(dist_mat(:,2),100*stats.MUSIC_DNN_STD,'x')
hold off
grid on
grid minor
title('Standard Deviation','interpreter','latex')
xlabel('Multi-Path Distance (m)','interpreter','latex')
ylabel('Standard Deviation (cm)','interpreter','latex')
legend('MUSIC',...
    'MUSIC-DNN',...
    'interpreter','latex',...
    'location','best');

% Mean Error
subplot(2,2,3);
hold on
scatter(dist_mat(:,2),abs(100*stats.MUSIC_mean_error),'o')
scatter(dist_mat(:,2),abs(100*stats.MUSIC_DNN_mean_error),'x')
scatter(dist_mat(:,2),abs(0*stats.MUSIC_DNN_mean_error),'k.')
hold off
grid on
grid minor
title('Absolute Distance Estimation Error','interpreter','latex')
xlabel('Multi-Path Distance (m)','interpreter','latex')
ylabel('Estimation Error (cm)','interpreter','latex')
legend('MUSIC',...
    'MUSIC-DNN',...
    'Theoretical',...
    'interpreter','latex',...
    'location','best');

% RMSE
subplot(2,2,4);
hold on
scatter(dist_mat(:,2),100*stats.MUSIC_RMSE,'o')
scatter(dist_mat(:,2),100*stats.MUSIC_DNN_RMSE,'x')
hold off
grid on
grid minor
title('RMSE','interpreter','latex')
xlabel('Multi-Path Distance (m)','interpreter','latex')
ylabel('RMSE (cm)','interpreter','latex')
legend('MUSIC',...
    'MUSIC-DNN',...
    'interpreter','latex',...
    'location','best');

sgtitle(plot_filename,'interpreter','latex')

% Figure Sizing
p_width = 1280;
p_height = 550;
set(1,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])

figure(3);
clf(3)
hold on
cdfplot(100*abs(results.MUSIC_err));
cdfplot(100*abs(results.MUSIC_DNN_err));
hold off
title('CDF Comparisons','interpreter','latex')
xlabel('Absolute Ranging Error (cm)','interpreter','latex')
ylabel('CDF','interpreter','latex')
legend('MUSIC',...
    'MUSIC-DNN',...
    'interpreter','latex','location','best')

p_width = 600;
p_height = 300;
set(3,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
end

function struct_out = cont_struct(struct_in)
for i = 1:size(struct_in,2)
    % Make sure to do initialization
    if i == 1
        struct_vrs = fieldnames(struct_in); % Fieldnames
        struct_out = struct_in(1); % First set of variable values
    else
        for k = 1:length(struct_vrs)
            struct_out.(struct_vrs{k}) = [struct_out.(struct_vrs{k});struct_in(i).(struct_vrs{k})];
        end
    end
end
end