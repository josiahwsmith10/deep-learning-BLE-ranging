% This script will plot the results using the different parameters from
% Table 1.

clear
close all

addpath(genpath('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\export_fig Yair Altman\'))

%% Load in the cases

load('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\case1')
stats1 = stats;
cal_dist_theo1 = cal_dist_theo;
load('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\case2')
stats2 = stats;
cal_dist_theo2 = cal_dist_theo;
load('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\case3')
stats3 = stats;
cal_dist_theo3 = cal_dist_theo;
load('C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\data\case4')
stats4 = stats;
cal_dist_theo4 = cal_dist_theo;

clear stats cal_dist_theo

%% Plot

try
    close 1
catch
end

fig1 = figure(1);
% Case 1 Plot
subplot(2,2,1)
hold on
plot(1e6*intra_delay,100*stats1.LS_mean_error,'o')
plot(1e6*intra_delay,100*cal_dist_theo1,'k')
hold off
grid on
% grid minor
title('Case 1','interpreter','latex')
xlabel('$T_o$ ($\mu s$)','interpreter','latex')
ylabel('Error (cm)','interpreter','latex')
legend('Sim.',...
    'Theo.',...
    'interpreter','latex','location','best')

% Case 2 Plot
subplot(2,2,2)
hold on
plot(1e6*inter_delay,100*stats2.LS_mean_error,'o')
plot(1e6*inter_delay,100*cal_dist_theo2,'k')
hold off
grid on
% grid minor
title('Case 2','interpreter','latex')
xlabel('$T_r$ ($\mu s$)','interpreter','latex')
ylabel('Error (cm)','interpreter','latex')
legend('Sim.',...
    'Theo.',...
    'interpreter','latex','location','best')

% Case 3 Plot
subplot(2,2,3)
hold on
plot(delta_v,100*stats3.LS_mean_error,'o')
plot(delta_v,100*cal_dist_theo3,'k')
hold off
grid on
% grid minor
title('Case 3','interpreter','latex')
xlabel('$\Delta_v$ (m/s)','interpreter','latex')
ylabel('Error (cm)','interpreter','latex')
legend('Sim.',...
    'Theo.',...
    'interpreter','latex','location','best')

% Case 4 Plot
subplot(2,2,4)
hold on
plot(1e6*crystal_offset,100*stats4.LS_mean_error,'o')
plot(1e6*crystal_offset,100*cal_dist_theo4,'k')
hold off
grid on
% grid minor
title('Case 4','interpreter','latex')
xlabel('$\epsilon$ (ppm)','interpreter','latex')
ylabel('Error (cm)','interpreter','latex')
legend('Sim.',...
    'Theo.',...
    'interpreter','latex','location','best')

% p_width = 400;
% p_height = 400;
p_width = 600;
p_height = 400;
set(1,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])

export_fig(fig1,'C:\Users\jayso\Box\Classes\Broadband Digital Communications\Project\MATLAB\Data\mob_crystal_offset_plot',...
    '-eps','-c[10,nan,nan,nan]','-transparent')