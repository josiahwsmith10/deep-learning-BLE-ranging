% This script simulates the BLE MCPD signal and proposed algorithms
function [results,stats] =...
    proj_mcpd_sim(num_runs,dist_vec,att_vec,...
    SNR,dnn_network,MUSIC_S,rayleigh)
%% Inputs
dist_LOS = dist_vec(1); % Get the LOS distance from the input distances vector
c = 3e8;                % Propagation speed
delta_f = 1e6;          % Measurement spacing in frequency
n = -37:37;             % Non-decimated set of subcarrier indexes
Nptot = length(n);      % Number of channels
cal_dist = 0;

%% MUSIC Inputs
M = size(MUSIC_S,1); % Subarray size
nsig = 5;
taus = cal_dist/c:0.001/c:50/c; % Possible delays considered.
dists = taus * c;
threshold_MUSIC = -13; % MUISC significant peak threshold in dB

%% Simulate the expected channels and combine them.
CFR_data_array = gen_ch(num_runs,SNR,dist_vec,att_vec,...
    n,delta_f,c,rayleigh);

%% Do multiple trials
dist_MUSIC_array = zeros(num_runs,nsig);
dist_MUSIC_DNN_array = zeros(num_runs,nsig);
for k = 1:num_runs
    %% Get the CFR vector for the run we are looking at
    CFR_data = CFR_data_array(:,k);
    
    %% MUSIC-based Distance Estimation
    [dist_MUSIC,spec_dB,~] = music_fb_alg(M,nsig,Nptot,dists,cal_dist,threshold_MUSIC,CFR_data,MUSIC_S);
    
    %% MUSIC-DNN-based Distance Estimation
    [dist_MUSIC_DNN,spec_DNN_dB,~] = music_dnn_alg(M,nsig,Nptot,dists,cal_dist,threshold_MUSIC,CFR_data,dnn_network,MUSIC_S);
    
    %% Output the Data
    dist_MUSIC_array(k,:) = dist_MUSIC.' / 2;
    dist_MUSIC_DNN_array(k,:) = dist_MUSIC_DNN.' / 2;
end

%% Get the results
method = 'quartiles';

% MUSIC results
x = dist_MUSIC_array(:,1);
x = x(~isoutlier(x,method));
results.MUSIC = x;
results.MUSIC_err = x - dist_LOS;

% MUSIC_DNN results
x = dist_MUSIC_DNN_array(:,1);
x = x(~isoutlier(x,method));
results.MUSIC_DNN = x;
results.MUSIC_DNN_err = x - dist_LOS;

%% Get statistics of the results array
% MUSIC Statistics
x = dist_MUSIC_array(:,1);
x = x(~isoutlier(x,method));
[~,C] = kmeans(x,3);
stats.MUSIC_mean = trimmean(C,30);
stats.MUSIC_STD = std(x);
stats.MUSIC_mean_error = stats.MUSIC_mean - dist_LOS;
stats.MUSIC_RMSE = sqrt(mean(abs(x - dist_LOS).^2));

% MUSIC_DNN Statistics
x = dist_MUSIC_DNN_array(:,1);
x = x(~isoutlier(x,method));
[~,C] = kmeans(x,3);
stats.MUSIC_DNN_mean = trimmean(C,30);

stats.MUSIC_DNN_STD = std(x);
stats.MUSIC_DNN_mean_error = stats.MUSIC_DNN_mean - dist_LOS;
stats.MUSIC_DNN_RMSE = sqrt(mean(abs(x - dist_LOS).^2));
