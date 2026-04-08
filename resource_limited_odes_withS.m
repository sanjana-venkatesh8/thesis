global light_exp_times light_exp_lengths;
light_exp_times = 0:100:500;%[10, 30, 50, 70, 90];
light_exp_lengths = 50 * ones(1, length(light_exp_times));

function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########  
    k_delay = 0.005;
    r_max = 0.05; r_min = 1e-4;
    P_max = 10e7;

    g_pdc1 = 3; r_cons = 2; R_p = 8;
    g_s = 1; d_s = 0.5; R_s = 10;
    g_N = 0.18; d_N = 0.05;
    g_c = 0.08;
    k_m_off = 0.01;
    %################################

    pdc1 = x(1);
    s = x(2);
    N = x(3);
    C = x(4);
    mem = x(5);
    rm = x(6);

    cell_logistic = 1;%(1 - (N + C)/P_max);

    dNdt = ((g_N) / (g_N*s + 1) * pdc1 - d_N * N) * cell_logistic;
    dsdt = (g_s - d_s) * N * (1 - s/(N * R_s));
    d_pdc1_dt = ((light(t)*g_pdc1) / (light(t)*g_pdc1*s + 1) * N - r_cons * dNdt) * (1 - pdc1/(N * R_p));
    dCdt = ((rm) / (rm * s + 1) * N - d_N * C + g_c * C) * cell_logistic;
    d_mem_dt = (1-light(t)) - k_m_off*mem;
    d_rm_dt = (k_delay*mem)*(r_max*(1-light(t)) - rm + r_min); % cannot use this to show slowed mutation rate

    dxdt = [d_pdc1_dt; 0; dNdt; dCdt; d_mem_dt; d_rm_dt];
    % dxdt = [d_pdc1_dt; dsdt; dNdt; dCdt; d_mem_dt; d_rm_dt];
end

function y = light(t)
    %########## PARAMETERS ##########
    global light_exp_times light_exp_lengths;
    %################################

    y = 0;
    for i = 1:length(light_exp_times)
        y = y + heaviside(t - light_exp_times(i)) ...
            - heaviside(t- (light_exp_times(i) + light_exp_lengths(i)));
    end
end

%########## INITIAL CONDITIONS ##########
pdc1_0 = 0.15; s_0 = 0;
N_0 = 1200; C_0 = 0;
mem0 = 0; rm0 = 0.02;
%########################################

[t,y] = ode45(@population_odes,[0 500],[pdc1_0; s_0; N_0; C_0; mem0; rm0]);

titles = ["PDC concentration vs time", "Chromatin binding protein concentration vs. time", "N cell population vs time", "Cheater population vs. time", "Total time in darkness vs. time", "Mutation rate vs. time"];
figure(1); hold on;
for i = 1:size(y, 2)
    subplot(size(y, 2), 1, i); hold on;
    plot(t, y(:, i), LineWidth=2); title(titles(i));
    xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
end

figure(2); hold on; plot(t, y(:, 3) + y(:, 4), LineWidth=2); title("total cell population vs. time")
plot(t, y(:, 3), LineWidth=2, LineStyle="--");
plot(t, y(:, 4), LineWidth=2, LineStyle='--');
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

%% PLOT EXP GROWTH TIME CONSTANT VS TIME (ln(N2) - ln(N1))/(t2 - t1)
N = y(:, 3) + y(:, 4);
mu = diff(log(N)) ./ diff(t);
t_mid = (t(1:end-1) + t(2:end)) / 2; % midpoints for plotting

figure;
plot(t_mid, mu, LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
xlabel('Time (h)');
ylabel('\mu (h^{-1})');
title('Population growth rate vs time');

% t_mid_2 = (t_mid(1:end-1) + t_mid(2:end)) / 2;
% figure; plot(t_mid_2, diff(mu) ./ diff(t_mid), LineWidth=2);
% xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

% %% plot ratio of pdc1 to cell count vs. time
% figure;plot(t, y(:,1) ./ y(:, 2), LineWidth=2); title("PDC1:cell count vs. time")