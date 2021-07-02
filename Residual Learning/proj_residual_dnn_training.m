%% Set Paramters
sim.att_min = 0.01;
sim.att_max = sqrt(10^(-3/10));

sim.crystal_offset_min = 0;
sim.crystal_offset_max = 0;

sim.SNR_min = 13;
sim.SNR_max = 17;

sim.intra_delay = 500e-6;
sim.Nptot = 75;

% MUSIC Parameters
sim.music.M = 30;
sim.music.c = physconst('lightspeed');
sim.music.nsig = 13;
sim.music.dists = 0:1e-3:50;
sim.music.threshold = -13;

% Generation Parameters
sim.numSamples = 16384;

sim.isRayleigh = true;

%% Generate Data
sim = generate_data(sim);

%% Examine the Data
examine_data(sim);

%% Create Network
net = initialize_net(sim);

%% Setup Training Data Features (X) and Labels (Y)
net.X = cat(3,...
    permute(real(sim.noisy_CFR),[1,3,4,2]),...
    permute(imag(sim.noisy_CFR),[1,3,4,2]));
net.Y = cat(3,...
    permute(real(sim.noise),[1,3,4,2]),...
    permute(imag(sim.noise),[1,3,4,2]));

%% Train Network
net.net = trainNetwork(net.X,net.Y,net.layers,net.options);

%% Test Network
sim_test = sim;
sim_test.numSamples = 1;
sim_test.music.nsig = 13;
sim_test = generate_data(sim_test);
sim_test.predX = cat(3,real(sim_test.noisy_CFR),imag(sim_test.noisy_CFR));
sim_test.pred = net.net.predict(sim_test.predX);
sim_test.pred = sim_test.pred(:,:,1) + 1j*sim_test.pred(:,:,2);
[sim_test.pred_music,sim_test.pred_pw,~] = MUSIC_pspec(sim_test,sim_test.music.M,sim_test.music.nsig,sim_test.Nptot,sim_test.noisy_CFR - sim_test.pred);
[sim_test.noisy_music,sim_test.noisy_pw,~] = MUSIC_pspec(sim_test,sim_test.music.M,sim_test.music.nsig,sim_test.Nptot,sim_test.noisy_CFR);
[sim_test.ideal_music,sim_test.ideal_pw,~] = MUSIC_pspec(sim_test,sim_test.music.M,sim_test.music.nsig,sim_test.Nptot,sim_test.ideal_CFR);
sim_test.lim = [min(sim_test.noisy_music)-3,max(sim_test.ideal_music)+3];

test_network(sim_test)

%% Save Network and MUSIC Steering Vectors
save_net(sim,net)

%% Functions
function examine_data(sim)
sim.indRand = randi(sim.numSamples);

sim.sgtitle = "Simulated: SNR = "+sim.SNR(sim.indRand)+" dB, Dists = "+num2str(sim.dist{sim.indRand})+" m, Path Magnitudes: "+num2str(sim.att{sim.indRand});

figure(12)

subplot(221)
plot(sim.music.dists,MUSIC_pspec(sim,sim.music.M,sim.music.nsig,sim.Nptot,sim.noisy_CFR(:,sim.indRand)));
title("MUSIC Spectrum of Noisy CFR",'interpreter','latex')
xlabel("Distance (m)",'interpreter','latex')
ylabel("Magnitude (dB)",'interpreter','latex')
grid on
grid minor

subplot(222)
plot(sim.music.dists,MUSIC_pspec(sim,sim.music.M,sim.music.nsig,sim.Nptot,sim.ideal_CFR(:,sim.indRand)));
title("MUSIC Spectrum of Noiseless CFR with " + sim.nsig(sim.indRand) + " targets",'interpreter','latex')
xlabel("Distance (m)",'interpreter','latex')
ylabel("Magnitude (dB)",'interpreter','latex')
grid on
grid minor

subplot(223)
plot(real(sim.noisy_CFR(:,sim.indRand)))
hold on
plot(real(sim.ideal_CFR(:,sim.indRand)))
plot(real(sim.noisy_CFR(:,sim.indRand) - sim.noise(:,sim.indRand)),'--')
hold off
xlim([1,sim.Nptot])
legend("Noisy","Ideal","Noise-Subtracted",'interpreter','latex')
title("Channel CFR Comparison",'interpreter','latex')
xlabel("Subchannel Index",'interpreter','latex')
ylabel("Magnitude",'interpreter','latex')
grid on
grid minor

subplot(224)
plot(sim.music.dists,MUSIC_pspec(sim,sim.music.M,sim.music.nsig,sim.Nptot,sim.noisy_CFR(:,sim.indRand) - sim.noise(:,sim.indRand)));
title("MUSIC Spectrum of Noise-Subtracted Ground-Truth CFR",'interpreter','latex')
xlabel("Distance (m)",'interpreter','latex')
ylabel("Magnitude (dB)",'interpreter','latex')
grid on
grid minor

sgtitle(sim.sgtitle,'interpreter','latex')

sgtitle(sim.sgtitle,'interpreter','latex')
p_width = 1280;
p_height = 550;
set(12,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
end

function net = initialize_net(sim)
net.layers = [ ...
    imageInputLayer([size(sim.noisy_CFR,1),1,2],"Normalization","none")
    
    convolution2dLayer([16,1],256,"Padding","same")
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer([13,1],256,"Padding","same")
    batchNormalizationLayer
    tanhLayer
    
    convolution2dLayer([8,1],128,"Padding","same")
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer([4,1],64,"Padding","same")
    batchNormalizationLayer
    tanhLayer
    
    convolution2dLayer([4,1],32,"Padding","same")
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer([4,1],32,"Padding","same")
    batchNormalizationLayer
    tanhLayer
    
    convolution2dLayer([4,1],32,"Padding","same")
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer([4,1],16,"Padding","same")
    batchNormalizationLayer
    tanhLayer
    
    convolution2dLayer([16,1],2,"Padding","same")
    
    regressionLayer
    ];

net.options = trainingOptions("adam", ...
    "MaxEpochs",1000, ...
    "InitialLearnRate",0.1,...
    "MiniBatchSize",64, ...
    "Plots","training-progress", ...
    "Verbose",true, ...
    "LearnRateSchedule","piecewise", ...
    "LearnRateDropFactor",0.9, ...
    "LearnRateDropPeriod",3, ...
    "ExecutionEnvironment","gpu", ...
    "Shuffle","every-epoch");
end

function sim = generate_data(sim)
delta_f = 1e6;
n = (-(sim.Nptot-1)/2:(sim.Nptot-1)/2).';
f_k = delta_f * n;
c = sim.music.c;

% Initialize holding variables
sim.ideal_CFR = zeros(sim.Nptot,sim.numSamples);
sim.noisy_CFR = zeros(sim.Nptot,sim.numSamples);
sim.noise = zeros(sim.Nptot,sim.numSamples);

% Precompute values
sim.nsig = 2*ones(1,sim.numSamples); % number of signals for each sample
sim.dist = cell(1,sim.numSamples);
sim.att = cell(1,sim.numSamples);
sim.eta_i = zeros(1,sim.numSamples);
sim.eta_r = zeros(1,sim.numSamples);
sim.SNR = zeros(1,sim.numSamples);
sim.f_k_i = zeros(sim.Nptot,sim.numSamples);
sim.f_k_r = zeros(sim.Nptot,sim.numSamples);
sim.T_k_i = zeros(1,sim.numSamples);
sim.T_k_r = zeros(1,sim.numSamples);
sim.time_offset = zeros(1,sim.numSamples);
sim.theta = zeros(sim.Nptot,sim.numSamples);

big_step = linspace(2,12,sqrt(sim.numSamples));
small_step = linspace(0,10,sqrt(sim.numSamples));
for indBig = 1:sqrt(sim.numSamples)
    for indSmall = 1:sqrt(sim.numSamples)
        sim.dist{(indBig-1)*sqrt(sim.numSamples)+indSmall} = [big_step(indBig),big_step(indBig)+small_step(indSmall)];
    end
end

for indSample = 1:sim.numSamples
    if sim.isRayleigh
        sim.att{indSample} = randn(1,sim.nsig(indSample)) + 1j*randn(1,sim.nsig(indSample));
    else
        sim.att{indSample} = [1,sim.att_min + (sim.att_max - sim.att_min)*rand(1,sim.nsig(indSample)-1)];
        sim.att{indSample} = sort(sim.att{indSample},'descend');
    end
    
    sim.eta_i(indSample) = sim.crystal_offset_min + (sim.crystal_offset_max - sim.crystal_offset_min)*rand();
    sim.eta_r(indSample) = sim.crystal_offset_min + (sim.crystal_offset_max - sim.crystal_offset_min)*rand();
    sim.SNR(indSample) = sim.SNR_min + (sim.SNR_max - sim.SNR_min)*rand();
    sim.f_k_i(:,indSample) = (sim.eta_i(indSample) + 1) * f_k;
    sim.f_k_r(:,indSample) = (sim.eta_r(indSample) + 1) * f_k;
    sim.T_k_i(indSample) = (sim.eta_i(indSample) + 1) * sim.intra_delay;
    sim.T_k_r(indSample) = (sim.eta_r(indSample) + 1) * 0;
    sim.time_offset(indSample) = sim.T_k_i(indSample) - sim.T_k_r(indSample);
    sim.theta(:,indSample) = pi*randn(sim.Nptot,1);
end

% MUSIC precompute steering vectors
sim.music.S = zeros(sim.music.M,length(sim.music.dists));
for indDist = 1:length(sim.music.dists)
    sim.music.S(:,indDist) = exp(-1j*2*pi*delta_f*sim.music.dists(indDist)/c*n(1:sim.music.M));
end

% Simulate channels
for indSample = 1:sim.numSamples
    att_vec = sim.att{indSample};
    dist_vec = sim.dist{indSample};
    time_offset = sim.time_offset(indSample);
    theta = sim.theta(indSample);
    f_k_i = sim.f_k_i(:,indSample);
    f_k_r = sim.f_k_r(:,indSample);
    SNR = sim.SNR(indSample);
    
    % Compute the channel responses
    channel_i = sum(att_vec .* exp(-1j*(2*pi*f_k_i.*(dist_vec/c - time_offset) - theta)),2);
    channel_r = sum(att_vec .* exp(-1j*(2*pi*f_k_r.*(dist_vec/c + time_offset) + theta)),2);
    
    % Get average signal energy for noise calculation
    Eave = mean([(1/sim.Nptot)*sum(abs(channel_i).^2),(1/sim.Nptot)*sum(abs(channel_r).^2)]);
    
    % Add noise for the initiator path
    noise_i = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(sim.Nptot,1) + 1i*randn(sim.Nptot,1));
    channel_i_n = channel_i + noise_i;
    
    % Add noise for the reflector path
    noise_r = sqrt(10^(-SNR/10)*(1/2)*Eave)*(randn(sim.Nptot,1) + 1i*randn(sim.Nptot,1));
    channel_r_n = channel_r + noise_r;
    
    % Set to sim struct
    sim.ideal_CFR(:,indSample) = channel_i .* channel_r;
    sim.noisy_CFR(:,indSample) = channel_i_n .* channel_r_n;
    
    % Noise remainder
    sim.noise(:,indSample) = channel_i .* noise_r + channel_r .* noise_i + noise_r .* noise_i;
end
end

function test_network(sim)
sim.sgtitle = "Simulated: SNR = "+sim.SNR+" dB, Dists = "+num2str(sim.dist{1})+" m, Path Magnitudes: "+num2str(sim.att{1});

figure(11)
clf
subplot(221)
plot(sim.music.dists,sim.noisy_music)
hold on
plot(sim.music.dists,sim.pred_music)
hold off
ylim(sim.lim)
title("MUSIC Spectrum of Noisy/Noise-Subtrated CFRs",'interpreter','latex')
legend("Noisy","Noise-Subtracted",'interpreter','latex')
xlim([sim.music.dists(1),sim.music.dists(end)])
xlabel("Distance (m)",'interpreter','latex')
ylabel("Magnitude (dB)",'interpreter','latex')
grid on
grid minor

subplot(222)
plot(sim.music.dists,sim.ideal_music)
ylim(sim.lim)
title("MUSIC Spectrum of Noiseless CFR",'interpreter','latex')
xlim([sim.music.dists(1),sim.music.dists(end)])
xlabel("Distance (m)",'interpreter','latex')
ylabel("Magnitude (dB)",'interpreter','latex')
grid on
grid minor

subplot(223)
plot(real(sim.noise))
hold on
plot(real(sim.pred))
hold off
xlim([1,sim.Nptot])
legend("Actual Noise","Predicted Noise",'interpreter','latex')
title("Ground Truth vs. Predicted",'interpreter','latex')
xlabel("Subchannel Index",'interpreter','latex')
ylabel("Magnitude",'interpreter','latex')
grid on
grid minor

subplot(224)
plot(real(sim.ideal_CFR))
hold on
plot(real(sim.noisy_CFR - sim.pred))
plot(real(sim.noisy_CFR))
hold off
xlim([1,sim.Nptot])
legend("Ideal Signal","Noise-Subtracted","Noisy",'interpreter','latex')
title("Ground Truth vs. Predicted",'interpreter','latex')
xlabel("Subchannel Index",'interpreter','latex')
ylabel("Magnitude",'interpreter','latex')
grid on
grid minor

sgtitle(sim.sgtitle,'interpreter','latex')
p_width = 1280;
p_height = 550;
set(11,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
end

function [pspec,peak_width,dist_MUSIC] = MUSIC_pspec(sim,M,nsig,Nptot,CFR_data)
% MATLAB's Smoothing Method for square covariance matrix. Need square
% covariance matrix for forward-backward averaging.
L = Nptot - M + 1;
RSM = complex(zeros(M,M));
for i = 1:L
    subArray = CFR_data(i:i+M-1);
    RSM = subArray*subArray' + RSM;
end
RSM = RSM / L;

J = fliplr(eye(M));
R = (RSM+J*conj(RSM)*J)/2; % Forward-backward averaging

% Compute the eigenvectors and eigenvalues
[eigenvects,eigenvals_diag,~] = svd(R);
eigenvals = diag(real(eigenvals_diag)); % column vector of eigenvalues

% Sort eigenvectors
[~,index] = sort(eigenvals,'descend');
eigenvects = eigenvects(:,index);

% Get the signal eigenvectors
% sig_eigenvects = eigenvects(:,1:nsig_mdl);

% Get the noise eigenvectors
noise_eigenvects = eigenvects(:,nsig+1:end);

% Compute MUSIC spectrum
spec_dens = abs((sim.music.S'*noise_eigenvects)).^2;
spec_den = sum(spec_dens,2);
pspec = (1./spec_den).';
pspec = 10*log10(pspec);

if nargout > 1
    spec_T = max(pspec + sim.music.threshold);
    % Find Delays from the peaks of the spectrum
    [~,locs,peak_width,~] = findpeaks(pspec(1,:),'SortStr','descend','MinPeakHeight',spec_T,...
        'WidthReference','halfprom'); % Sorted by peakage magnitude. Excluding negative delays.
    
    % Handling if no peaks found
    if isempty(locs)
        
        [~,locs] = max(pspec(1,:)); % Use max if no peaks found
        peak_width = inf; % Set peakwidth to inf
        
    else
        
        locs_check = min(sim.music.nsig,length(locs)); % Incase we don't detect at least nsig peaks
        [locs,sort_ind] = sort(locs(1:locs_check),'ascend'); % Sorted by smallest to largest delay.
        
        % Format the vector of peak widths
        peak_width = peak_width(sort_ind);
        peak_width = peak_width(1);
    end
    dist_MUSIC = sim.music.dists(locs)/2;
end
end

function save_net(sim,net)
S = sim.music.S;
if sim.isRayleigh
    save net_AWGN_Rayleigh net S
else
    save net_AWGN net S
end
end
