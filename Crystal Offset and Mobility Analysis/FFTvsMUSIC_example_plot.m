% This script will combine the example plots for SNR = 15 dB and SNR = 25
% dB.

clear
close all

%% Load Saved figures

fig15 = hgload('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\FFTvsMUSIC_example_SNR15_2m6.5m_15dB.fig');
fig25 = hgload('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\FFTvsMUSIC_example_SNR15_2m6.5m_25dB.fig');

%% Prepare subplots
fig3 = figure;
h(1) = subplot(2,1,1);
h(2) = subplot(2,1,2);

%% Paste figures on the subplots

copyobj(allchild(get(fig15,'CurrentAxes')),h(1));
copyobj(allchild(get(fig25,'CurrentAxes')),h(2));

%% Add plot parameters

fig3;
subplot(2,1,1);
grid on
% grid minor
title({'Two-Way FFT vs. MUSIC Spectrum Comparisons','SNR = 15 dB'},'interpreter','latex')
xlabel('Distance (m)','interpreter','latex')
ylabel('Magnitude (dB)','interpreter','latex')
legend('MUSIC',...
    'FFT',...
    'True Multi-Path Distances',...
    'interpreter','latex','location','best')

subplot(2,1,2);
grid on
% grid minor
title('SNR = 25 dB','interpreter','latex')
xlabel('Distance (m)','interpreter','latex')
ylabel('Magnitude (dB)','interpreter','latex')
legend('MUSIC',...
    'FFT',...
    'True Multi-Path Distances',...
    'interpreter','latex','location','best')

%% Resize the figure

% p_width = 600;
% p_height = 520;
p_width = 520;
p_height = 520;
set(fig3,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])


%% Explort with export_fig

addpath(genpath('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\export_fig Yair Altman\'))

export_fig(fig3,'C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\Data\FFTvsMUSIC_example_latex',...
    '-eps','-c[0,nan,nan,nan]','-transparent')