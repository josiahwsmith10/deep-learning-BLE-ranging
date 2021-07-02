% Quantifying errors

%%
c = 3e8;
eta_r = -10e-6;
eta_i = 10e-6;
f_0 = 2.404e9;
delta_f = 1e6;
K_f = 75;

alpha_f = -(c/2) * (eta_i - eta_r)^2 * (f_0/delta_f + K_f);
alpha_o = -(c/2) * (1+eta_i)*(eta_i - eta_r);

% alpha_o / alpha_f

%% Doppler

f_0 = 2.404e9;
f_Kf = 2.478e9;
v_i = 3;
v_r = -3;

f_k_i = (1 + eta_i)*f_Kf;
f_k_r = (1 + eta_r)*f_Kf;

f_k_i_d = (f_k_r / c)*(v_i - v_r);
f_k_r_d = (f_k_i / c)*(v_r - v_i);

a_f = (c / (2*delta_f))*(f_k_i_d - f_k_r_d)*(eta_r - eta_i);
a_o = -(1/2)*((1+eta_i)^2 + (1+eta_i)*(1+eta_r))*(v_i - v_r);

alpha_f_new = alpha_f + a_f;
alpha_o_new = alpha_o + a_o;

alpha_f_change = (alpha_f_new - alpha_f) / alpha_f
alpha_o_change = (alpha_o_new - alpha_o) / alpha_o