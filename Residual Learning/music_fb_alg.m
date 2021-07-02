% MUSIC distance estimation algorithm with forwward-backward averaging and
% MATLAB's covariance matrix generation which results in a square matrix.

function [dist_MUSIC,spec_dB,peakwidth] = music_fb_alg(M,nsig,Nptot,dists,cal_dist,...
    threshold,CFR_data,S)
num_ant = size(CFR_data,2);

if num_ant > 1
    % Combining across antenna measurements
    L = Nptot - M + 1;
    num_ant = size(CFR_data,2);
    Rn = complex(zeros(Nptot,Nptot));
    for k = 1:num_ant
        Rn = CFR_data(:,k)*CFR_data(:,k)' + Rn;
    end
    Rn = Rn / num_ant;
    
    % Smoothing
    RSM = complex(zeros(M,M)); % Correlation matrix
    for k = 1:L
        RSM = RSM + Rn(k:k+M-1,k:k+M-1);
    end
    RSM = RSM / L;
else
    % MATLAB's Smoothing Method for square covariance matrix. Need square
    % covariance matrix for forward-backward averaging.
    L = Nptot - M + 1;
    RSM = complex(zeros(M,M));
    for i = 1:L
        subArray = CFR_data(i:i+M-1);
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
% spec_T = MUSIC_peaks(1) * 0.05;
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
dist_MUSIC = dists(locs).' - cal_dist;

% Include this so that the plot title doesn't give an error if we only
% calculate one path.
if length(dist_MUSIC) < nsig && nsig > 1
    dist_MUSIC = [dist_MUSIC;nan(nsig-length(dist_MUSIC),1)];
end