% This function will generate either a two-way channel response in either
% AWGN only or Rayleigh fading with AWGN.

function CFR_data_array = gen_ch(num_runs,SNR,dist_vec,att_vec,...
    n,delta_f,c,rayleigh)
% Using a random attenuation vector for Rayleigh fading
if rayleigh == true
    clear att_vec
end

npaths = length(dist_vec);
Nptot = length(n);

CFR_data_array = zeros(Nptot,num_runs);

for k = 1:num_runs
    % Get random attenuations for Rayleigh
    if rayleigh == true
        % att_vec = abs(randn(1,npaths));
        att_vec = randn(1,npaths) + 1i*randn(1,npaths);
    end
    
    % Generate the initiator->reflector and reflector->initiator channel responses
    channel_comb_i = zeros(Nptot,1);
    channel_comb_r = zeros(Nptot,1);
    channels_i = zeros(Nptot,npaths);
    channels_r = zeros(Nptot,npaths);
    for i = 1:npaths
        channels_i(:,i) = att_vec(i)*exp(-1i*(2*pi*delta_f*(dist_vec(i)/(c))*n.'));
        channels_r(:,i) = att_vec(i)*exp(-1i*(2*pi*delta_f*(dist_vec(i)/(c))*n.'));
        
        channel_comb_i = channel_comb_i + channels_i(:,i);
        channel_comb_r = channel_comb_r + channels_r(:,i);
    end
    
    % Get average signal energy for noise calculation
    Eave = mean([(1/Nptot)*sum(abs(channel_comb_i).^2),(1/Nptot)*sum(abs(channel_comb_r).^2)]);
    
    % Add noise for the initiator -> reflector path
    noise(:,1) = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(Nptot,1) + 1i*randn(Nptot,1));
    channel_comb_rx(:,1) = channel_comb_i + noise(:,1);
    
    % Add noise for the reflector -> initiator path
    noise(:,2) = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(Nptot,1) + 1i*randn(Nptot,1));
    channel_comb_rx(:,2) = channel_comb_r + noise(:,2);
    
    % Multiply the two signals and make column vector
    CFR_data_array(:,k) = (channel_comb_rx(:,1) .* channel_comb_rx(:,2));
end
