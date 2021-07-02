% Control script for proj_mcpd_sim. Change the desired input parameters and
% the two for loops to observe the effect of different impairments.

clear
close all

%% Inputs

% addpath('C:\Users\jayso\Box\Wireless Lab\BLE Localization\MATLAB\')

% dist_mat = [10*ones(160,1),linspace(10.05,18,160).']; 
% att_mat = repmat([1,sqrt(10^(-3/10))],size(dist_mat,1),1); 
dist_mat = [2,6.5]; % Trial distances in meters. Each row is a trial. Each column is a path distance.
att_mat = [1,1]; % Path attenuation matrix for the distances
intra_delay = 0; % Intra or Inter delay - Change which this is in proj_mcpd_sim.m!
crystal_offset = 0; % Crystal offset between initiator and reflector
delta_v = 0; % Doppler velocity (v_i - v_r) in m/s
SNR = 25; % SNR to simulate in dB
num_runs = 1000; % Sets the number of trials to do for each distance
showplot = 0; % Indicator to show plots in the ranging script

%% Processing one distance at a time

figure_handle = 10:10:10*length(intra_delay)*size(dist_mat,1);
ctr = 1;
for i = 1:length(intra_delay)

    for k = 1:size(dist_mat,1)
        
        disp(['Trial: ',num2str(ctr),'/',num2str(size(dist_mat,1)*length(intra_delay))]);
        
        % Run the ranging script
        [results_all(ctr),stats_all(ctr),cal_dist_theo(ctr)] =...
            proj_mcpd_sim(num_runs,dist_mat(k,:),att_mat(k,:),intra_delay(i),...
            crystal_offset,delta_v,SNR,figure_handle(ctr),showplot);
        
        ctr = ctr + 1;
        
    end
    
end

%% Concatenate structures
results = cont_struct(results_all);
stats = cont_struct(stats_all);

%% Define a plot name

plot_filename = strcat("Simulated: $\epsilon$ = ",num2str(1e6*crystal_offset),...
    " ppm, $\Delta_v$ = ",num2str(delta_v)," Hz, SNR = ",num2str(SNR)," dB");

%% Plot Stats
try
    close 1
    close 2
    close 3
catch
end

figure(1);
% Mean Error across the desired parameter - used to look at the effect of
% changing intra-delay, inter-delay, crystal offset, or mobility.
hold on
plot(1e6*intra_delay,100*stats.LS_mean_error,'o')
% plot(1e6*trial_delays,100*stats.MUSIC_mean_error,'+')
% plot(1e6*trial_delays,100*stats.FFT_mean_error,'x')
plot(1e6*intra_delay,100*cal_dist_theo,'color','k')
hold off
grid on
grid minor
title('Distance Estimation Error','interpreter','latex')
xlabel('Reflector Delay ($\mu s$)','interpreter','latex')
ylabel('Estimation Error (cm)','interpreter','latex')
legend('Estimated',...
    'Theoretical',...
    'interpreter','latex',...
    'location','best');

% Use to analyze the effect for different multipath distances.
if size(dist_mat,2) > 1
    figure(2);

    % Estimated distances
    subplot(2,2,1);
    hold on
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.MUSIC_mean)
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.FFT_mean)
    hold off
    grid on
    grid minor
    title({plot_filename,'Distance Estimations'},'interpreter','latex')
    xlabel('Multipath Distance $\Delta$ (cm)','interpreter','latex')
    ylabel('Estimated Distance (cm)','interpreter','latex')
    legend({'MUSIC',...
        'FFT',...
        'Theoretical'},...
        'interpreter','latex',...
        'location','best');

    % Standard Deviations
    subplot(2,2,2);
    hold on
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.MUSIC_STD)
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.FFT_STD)
    hold off
    grid on
    grid minor
    title('Standard Deviation','interpreter','latex')
    xlabel('Multipath Distance $\Delta$ (cm)','interpreter','latex')
    ylabel('Standard Deviation (cm)','interpreter','latex')
    legend('MUSIC',...
        'FFT',...
        'interpreter','latex',...
        'location','best');

    % Mean Error
    subplot(2,2,3);
    hold on
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.MUSIC_mean_error)
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.FFT_mean_error)
    hold off
    grid on
    grid minor
    title('Distance Estimation Error','interpreter','latex')
    xlabel('Multipath Distance $\Delta$ (cm)','interpreter','latex')
    ylabel('Estimation Error (cm)','interpreter','latex')
    legend('MUSIC',...
        'FFT',...
        'Theoretical',...
        'interpreter','latex',...
        'location','best');

    % RMSE
    subplot(2,2,4);
    hold on
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.MUSIC_RMSE)
    plot(100*(dist_mat(:,2)-dist_mat(1,1)),100*stats.FFT_RMSE)
    hold off
    grid on
    grid minor
    title('RMSE','interpreter','latex')
    xlabel('Multipath Distance $\Delta$ (cm)','interpreter','latex')
    ylabel('RMSE (cm)','interpreter','latex')
    legend('MUSIC',...
        'FFT',...
        'interpreter','latex',...
        'location','best');

    % Figure Sizing
    p_width = 1280;
    p_height = 550;
    set(2,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])


    figure(3);
    hold on
    cdfplot(100*abs(results.MUSIC_err));
    cdfplot(100*abs(results.FFT_err));
    title('CDF Comparisons','interpreter','latex')
    xlabel('Absolute Ranging Error (cm)','interpreter','latex')
    ylabel('CDF','interpreter','latex')
    legend('MUSIC',...
        'FFT',...
        'interpreter','latex','location','best')

    p_width = 1280;
    p_height = 550;
    set(3,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
end