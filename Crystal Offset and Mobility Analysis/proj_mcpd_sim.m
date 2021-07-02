% This script simulates the BLE MCPD signal and algorithms.

% clear
% close all

function [results,stats,cal_dist_theo] =...
        proj_mcpd_sim(num_runs,dist_vec,att_vec,delay_in,eta_i,delta_v,...
        SNR,figure_handle,showplot)

% figure_handles = [double(uint32(intra_delay*10000 + 10)),double(uint32(intra_delay*10000 + 20))]; 
    
%% Inputs

dist_LOS = dist_vec(1); % Get the LOS distance from the input distances vector
npaths = length(dist_vec); % Number of one-way paths simulated
c = 3e8; % Propagation speed
delta_f = 1e6; % Measurement spacing in frequency
f_0 = 2.404e9; % First channel in Hz
n = -37:37; % Non-decimated set of subcarrier indexes
Nptot = length(n);

% Clock offset for the reflector. Set as zero and vary the initiator clock
% offset for simplicity.
eta_r = 0; 

% Values from the BLE paper
% eta_i = 1 + 4e-6;
% eta_r = 1 + 7.25e-6;

% Include effect of crystal offset on carrier frequencies
f_k = (f_0 + (0:delta_f:(Nptot-1)*delta_f)).';
f_k_i = (eta_i + 1) * f_k; % Initiator frequency step with offset
f_k_r = (eta_r + 1) * f_k; % Reflector frequency step with offset

% Include effect of mobility
f_d_k_i = (((1 + eta_r)*f_k)/c) * (delta_v); % Doppler frequency shift at initiator
f_d_k_r = (((1 + eta_i)*f_k)/c) * (-delta_v); % Doppler frequency shift at reflector
f_k_i = f_k_i + f_d_k_i;
f_k_r = f_k_r + f_d_k_r;

T_r = 1e-3; % Inter-delay
T_o = delay_in; % Intra-delay
T_f = T_o + T_r; % Total delay for one measurement
t_k = linspace(0,(Nptot-1)*T_f,Nptot);
T_k_i = (eta_i + 1) * (t_k + T_o); % Received times at the initiator with offset
T_k_r = (eta_r + 1) * t_k; % Received times at the reflector

time_offset = (T_k_i - T_k_r).'; % Time offset (delta_t in 2019 BLE paper)

rng('shuffle')
theta = pi*randn(Nptot,1); % PLL offset per subcarrier
% theta = 0;

% Theoretical one-way calibration distance
alpha_f = -(c/(2*delta_f))*((eta_i - eta_r)^2 * (f_0 + (Nptot-1)*delta_f) - ...
    (f_d_k_i(end) - f_d_k_r(end))*(eta_r - eta_i));
alpha_o = -(c/2)*(1 + eta_i)*(eta_i - eta_r) - ...
    (1/2)*((1 + eta_i)^2 + (1 + eta_i)*(1 + eta_r))*delta_v;

cal_dist_theo = alpha_f*T_f + alpha_o*T_o;

% cal_dist_theo = -c/2 * (eta_i - eta_r) * T_o; 
cal_dist = 0;

%% MUSIC Inputs
M = 30; % Subarray size
nsig = 5;
taus = cal_dist/c:0.01/c:18/c; % Possible delays considered. 
dists = taus * c;
threshold_MUSIC = -13; % MUISC significant peak threshold in dB

%% FFT Inputs
% FFT_size = 2^23;
FFT_size = 2^13;
w = -pi:2*pi/FFT_size:pi - 2*pi/FFT_size;
FFT_dists = (w * c) ./ -(2*pi*delta_f);

% Get desired FFT indexes to search
FFT_indexes_1 = FFT_dists >= cal_dist;
FFT_indexes_2 = FFT_dists <= 18;
FFT_indexes = FFT_indexes_1 + FFT_indexes_2 == 2;
FFT_dists_indexed = FFT_dists(FFT_indexes);

threshold_FFT = -13; % FFT peak threshold in dB

%% Simulate the expected channels and combine them.

% Create the CFR array before processing the trials
CFR_data_array = complex(zeros(Nptot,num_runs));
for k = 1:num_runs
    
    % Generate the initiator->reflector and reflector->initiator channel responses
    channel_comb_i = zeros(Nptot,1);
    channel_comb_r = zeros(Nptot,1);
    channels_i = complex(zeros(Nptot,npaths));
    channels_r = complex(zeros(Nptot,npaths));
    for i = 1:npaths
        channels_i(:,i) = att_vec(i)*exp(-1i*(2*pi*f_k_i.*(dist_vec(i)/c - time_offset) - theta));
        channels_r(:,i) = att_vec(i)*exp(-1i*(2*pi*f_k_r.*(dist_vec(i)/c + time_offset) + theta));

        channel_comb_i = channel_comb_i + channels_i(:,i);
        channel_comb_r = channel_comb_r + channels_r(:,i);
    end
    
    % Get average signal energy for noise calculation
    Eave = mean([(1/Nptot)*sum(abs(channel_comb_i).^2),(1/Nptot)*sum(abs(channel_comb_r).^2)]);
    
    % Add noise for the initiator -> reflector path
    rng('shuffle')
    noise(:,1) = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(Nptot,1) + 1i*randn(Nptot,1));
    channel_comb_rx(:,1) = channel_comb_i + noise(:,1);

    % Add noise for the reflector -> initiator path
    rng('shuffle')
    noise(:,2) = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(Nptot,1) + 1i*randn(Nptot,1));
    channel_comb_rx(:,2) = channel_comb_r + noise(:,2);

    % Multiply the two signals and make column vector
    CFR_data_array(:,k) = (channel_comb_rx(:,1) .* channel_comb_rx(:,2));
end


%% Do multiple trials
dist_LS_array = zeros(num_runs,1);
dist_MUSIC_array = zeros(num_runs,nsig);
dist_FFT_array = zeros(num_runs,nsig);
for k = 1:num_runs

    %% Get the CFR vector for the run we are looking at
    CFR_data = CFR_data_array(:,k);
    
    %% Calculate the channel phase response for computing distance based on LS and plotting
    
    CPR_data = unwrap(angle(CFR_data));
   
    %% Lease Squares-based Distance estimation
    [dist_LS,~] = ls_alg(n,delta_f,cal_dist,c,CPR_data);
    
    %% MUSIC-based Distance Estimation
    [dist_MUSIC,spec_dB,~] = music_fb_alg(M,nsig,Nptot,delta_f,dists,c,n,cal_dist,threshold_MUSIC,CFR_data);
    
    %% FFT-based Distance Estimation
    [dist_FFT,FFT_spec_dB,~] = fft_alg(FFT_size,FFT_dists_indexed,FFT_indexes,nsig,...
    cal_dist,threshold_FFT,CFR_data);
    
    %% Both Magnitude and Phase Reponse Plots
    
    if showplot == 1
        figure(figure_handle);

        % Plot Magnitude Response
        subplot(2,2,1);
        hold on
        plot(n,10*log10(abs(CFR_data).^2),'.');
        hold off

        % Plot Phase Response
        subplot(2,2,3);
        hold on
        plot(n,CPR_data,'.');
        hold off

        % Plot MUSIC Spectrum
        subplot(2,2,2);
        hold on
        plot(dists-cal_dist,spec_dB)
        hold off

        % Plot FFT Spectrum
        subplot(2,2,4);
        hold on
        plot(FFT_dists_indexed-cal_dist,FFT_spec_dB(FFT_indexes))
        hold off
        
        
        % Compare MUSIC and FFT Plots
        figure(figure_handle + 1);
        hold on
        plot(dists-cal_dist,spec_dB)
        plot(FFT_dists_indexed-cal_dist,FFT_spec_dB(FFT_indexes))
        hold off
    end
    

    %% Output the Data
    
    dist_LS_array(k,:) = dist_LS / 2;
    dist_MUSIC_array(k,:) = dist_MUSIC.' / 2;
    dist_FFT_array(k,:) = dist_FFT.' / 2;

end

%% Get the results

% LS results
results.LS = dist_LS_array(:,1);
results.LS_err = dist_LS_array(:,1) - dist_LOS;

% MUSIC results
results.MUSIC = dist_MUSIC_array(:,1);
results.MUSIC_err = dist_MUSIC_array(:,1) - dist_LOS;

% FFT results
results.FFT = dist_FFT_array(:,1);
results.FFT_err = dist_FFT_array(:,1) - dist_LOS;

%% Get statistics of the results array

% LS Statistics
stats.LS_mean = mean(dist_LS_array(:,1));
stats.LS_STD = std(dist_LS_array(:,1));
stats.LS_mean_error = stats.LS_mean - dist_LOS;
stats.LS_RMSE = sqrt(mean(abs(dist_LS_array(:,1) - dist_LOS).^2));

% MUSIC Statistics
stats.MUSIC_mean = mean(dist_MUSIC_array(:,1));
stats.MUSIC_STD = std(dist_MUSIC_array(:,1));
stats.MUSIC_mean_error = stats.MUSIC_mean - dist_LOS;
stats.MUSIC_RMSE = sqrt(mean(abs(dist_MUSIC_array(:,1) - dist_LOS).^2));

% FFT Statistics
stats.FFT_mean = mean(dist_FFT_array(:,1));
stats.FFT_STD = std(dist_FFT_array(:,1));
stats.FFT_mean_error = stats.FFT_mean - dist_LOS;
stats.FFT_RMSE = sqrt(mean(abs(dist_FFT_array(:,1) - dist_LOS).^2));


%% Define the plot filename

plot_filename = strcat("Simulated: SNR = ",num2str(SNR)," dB, Dists = ",num2str(dist_vec)," m, Path Magnitudes: ",num2str(att_vec));


%% Both Magnitude and Phase Reponse Plots
if showplot == 1
    figure(figure_handle);
    % Plot Magnitude Response
    subplot(2,2,1);
    grid on
    grid minor
    % axis([-inf,inf,0,max(abs(CFR_data_array),[],'all')])
    title({plot_filename,'Two-Way Channel Magnitude Responses'},'interpreter','latex')
    xlabel('Subcahnnel Index', 'interpreter', 'latex')
    ylabel('Magnitude (dB)', 'interpreter', 'latex')
    % legend('Two-Way Magnitude Response Measurements',...
    % 'interpreter', 'latex', 'location', 'best')

    % Plot Phase Response
    subplot(2,2,3);
    grid on
    grid minor
    title("Two-Way Channel Phase Responses", 'interpreter', 'latex')
    xlabel('Subcahnnel Index', 'interpreter', 'latex')
    ylabel('Unwrapped Phase (rad.)', 'interpreter', 'latex')
    % legend('Two-Way Phase Response Measurements',...
    %     'Linear Fits',...
    % 'interpreter', 'latex', 'location', 'best')

    % Plot MUSIC Spectrum
    subplot(2,2,2);
    grid on
    grid minor
    title(['MUSIC, One-way Actual Dist: ',num2str(dist_LOS),...
        ' m, Calculated Dist: ',num2str(stats.MUSIC_mean)],'interpreter','latex')
    xlabel('Distance (m)','interpreter','latex')
    ylabel('Magnitude (dB)','interpreter','latex')

    % Plot FFT Spectrum
    subplot(2,2,4);
    grid on
    grid minor
    title(['FFT, One-way Actual Dist: ',num2str(dist_LOS),...
        ' m, Calculated Dist: ',num2str(stats.FFT_mean)],'interpreter','latex')
    xlabel('Distance (m)','interpreter','latex')
    ylabel('Magnitude (dB)','interpreter','latex')

    % Figure Sizing
    p_width = 1280;
    p_height = 550;
    set(figure_handle,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])

    figure(figure_handle + 1);
    % Plot FFT Spectrum
    xline(2*dist_vec(1),'color','k')
    xline((dist_vec(1) + dist_vec(2)),'color','k')
    xline(2*dist_vec(2),'color','k')
    grid on
    grid minor
    title(['Two-Way FFT vs. MUSIC Spectrum Comparison (SNR = ',num2str(SNR),' dB)'],'interpreter','latex')
    xlabel('Distance (m)','interpreter','latex')
    ylabel('Magnitude (dB)','interpreter','latex')
    legend('MUSIC',...
        'FFT',...
        'True Multi-Path Distances',...
        'interpreter','latex','location','best')

    % Figure Sizing
    p_width = 600;
    p_height = 500;
    set(figure_handle+1,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
end