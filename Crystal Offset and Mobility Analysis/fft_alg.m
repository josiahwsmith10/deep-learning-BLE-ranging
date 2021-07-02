% FFT Distance Estimation

function [dist_FFT,FFT_comp_spec_dB,peakprom] = fft_alg(FFT_size,FFT_dists_indexed,FFT_indexes,nsig,...
    cal_dist,threshold,CFR_data)
    
% Take FFT
FFT_comp_spec = abs(fftshift(fft(CFR_data,FFT_size))).^2;
FFT_comp_spec = FFT_comp_spec ./ max(FFT_comp_spec);

% Get spectrum in dB
FFT_comp_spec_dB = 10*log10(FFT_comp_spec);

% Define a threshold for FFT peak search
spec_T_FFT = max(FFT_comp_spec) * 10^(threshold/10);

% Get distances
[~,locs,~,peakprom] = findpeaks(FFT_comp_spec(FFT_indexes),'SortStr','descend',...
    'MinPeakHeight',spec_T_FFT);


if isempty(locs)
    
    [~,locs] = max(FFT_comp_spec(FFT_indexes)); % Use max if no peaks found
    
    dist_FFT = FFT_dists_indexed(locs) - cal_dist;
    
    peakprom = 0; % Set peakprom to 0
    
    warning('No FFT peak detected!')
    
else
    
    % In case we have less peaks than expected number of paths
    locs_check = min(nsig,length(locs)); 
    locs = locs(1:locs_check);
    
    dist_FFT = FFT_dists_indexed(locs) - cal_dist;
    [dist_FFT,sort_ind] = sort(dist_FFT,'ascend'); % Sort from smallest to largest distance
    dist_FFT = dist_FFT.'; % Make column vector
    
    % Format the vector of peak prominences
    peakprom = peakprom(sort_ind);
    peakprom = peakprom(1);

end

% Compensate if there are less peaks than number of paths
if length(dist_FFT) < nsig && nsig > 1
    dist_FFT = [dist_FFT;nan(nsig-length(dist_FFT),1)];
end