% Least Squares Distance Estimation by linear fit

function [dist_LS,CPR_fit] = ls_alg(SC_indexes,delta_f,cal_dist,c,CPR_data)

SC_indexes_col = SC_indexes.';
SC_indexes_col_LS = [SC_indexes.', ones(length(SC_indexes),1)];

% Only need compensated slope and linear fit
CPR_LS = inv(SC_indexes_col_LS.' * SC_indexes_col_LS) * SC_indexes_col_LS.' * CPR_data;
CPR_fit = (CPR_LS(1).*SC_indexes_col) + CPR_LS(2);

%% Calculate the delay using the combined phase slope
tau_LS = CPR_LS(1) ./ -(2*pi*delta_f);
dist_LS = (tau_LS .* c) - cal_dist;