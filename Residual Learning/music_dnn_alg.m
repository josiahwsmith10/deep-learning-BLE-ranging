% MUSIC distance estimation algorithm with forwward-backward averaging and
% MATLAB's covariance matrix generation which results in a square matrix.
function [dist_MUSIC_DNN,spec_dB,peakwidth] = music_dnn_alg(M,nsig,Nptot,dists,cal_dist,...
    threshold,CFR_data,dnn_network,S)
num_ant = size(CFR_data,2);

if num_ant > 1
    error("num_ant > 1 is not supported by DNN")
elseif Nptot ~= dnn_network.Layers(1).InputSize(1)
    error("Nptot must equal " + dnn_network.Layers(1).InputSize(1) + " for specified network")
else
    X = cat(3,real(CFR_data),imag(CFR_data));
    Y = dnn_network.predict(X);
    pred_noise = Y(:,:,1) + 1j*Y(:,:,2);
    CFR_data_denoised = CFR_data - pred_noise;
    
    % MATLAB's Smoothing Method for square covariance matrix. Need square
    % covariance matrix for forward-backward averaging.
    L = Nptot - M + 1;
    RSM = complex(zeros(M,M));
    for i = 1:L
        subArray = CFR_data_denoised(i:i+M-1);
        RSM = subArray*subArray' + RSM;
    end
    RSM = RSM / L;
end

J = fliplr(eye(M));
R = (RSM+J*conj(RSM)*J)/2; % Forward-backward averaging

% Compute the eigenvectors and eigenvalues
[eigenvects,eigenvals_diag,~] = svd(R);
eigenvals = diag(real(eigenvals_diag)); % column vector of eigenvalues

% Sort eigenvectors
[~,index] = sort(eigenvals,'descend');
eigenvects = eigenvects(:,index);

% Get the noise eigenvectors
noise_eigenvects = eigenvects(:,nsig+1:end);

% Compute MUSIC spectrum
spec_dens = abs((S'*noise_eigenvects)).^2;
spec_den = sum(spec_dens,2);
spec = (1./spec_den).';
spec_dB = 10*log10(spec);

% Define peak threshold
spec_T = max(spec(1,:)) * 10^(threshold/10);

% Find Delays from the peaks of the spectrum
[~,locs,peakwidth,~] = findpeaks(spec(1,:),'SortStr','descend','MinPeakHeight',spec_T,...
    'WidthReference','halfprom'); % Sorted by peakage magnitude. Excluding negative delays.

% Handling if no peaks found
if isempty(locs)
    [~,locs] = max(spec(1,:)); % Use max if no peaks found
    peakwidth = inf; % Set peakwidth to inf
else
    locs_check = min(nsig,length(locs)); % Incase we don't detect at least nsig peaks
    [locs,sort_ind] = sort(locs(1:locs_check),'ascend'); % Sorted by smallest to largest delay.
    
    % Format the vector of peak widths
    peakwidth = peakwidth(sort_ind);
    peakwidth = peakwidth(1);
end

% Calculate distance
dist_MUSIC_DNN = dists(locs).' - cal_dist;

% Include this so that the plot title doesn't give an error if we only
% calculate one path.
if length(dist_MUSIC_DNN) < nsig && nsig > 1
    dist_MUSIC_DNN = [dist_MUSIC_DNN;nan(nsig-length(dist_MUSIC_DNN),1)];
end
end