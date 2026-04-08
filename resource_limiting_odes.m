global light_exp_times light_exp_lengths;
light_exp_times = 0;%[10, 30, 50, 70, 90];
light_exp_lengths = 0 * ones(1, length(light_exp_times));

function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########    
    g_cells = 0.18; r_death = 0.05; cells_max = 1000000;
    g_pdc1 = 3;
    r_pdcCons = 2;
    R = 8;
    r_m = 5e-6; r_gm = 0.20;
    g_s = 2; R_s = 10;
    %################################

    pdc1 = x(1);
    cells = x(2);
    mem = x(3);
    M = x(4);
    S = x(5);

    dcellsdt = (g_cells*pdc1 * 1/S -r_death*cells)*(1-(cells+M)/cells_max);
    dpdc1dt = (light(t)*g_pdc1*cells * 1/S - dcellsdt * r_pdcCons) * (1 - pdc1/(cells*R));
    dmemdt = 1 - light(t); % total time in darkness
    dMdt = (r_m * mem * cells * 1/S + (r_gm - r_death) * M)*(1-(cells+M)/cells_max);
    dSdt = 0;%g_s * cells * (1 - S/(cells*R_s));

    dxdt = [dpdc1dt; dcellsdt; dmemdt; dMdt; dSdt];
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
pdc1_0 = 0.15; cells_0 = 1200; mem0 = 0; M_0 = 0; S_0 = 1;%1e-5;
%########################################

[t,y] = ode45(@population_odes,[0 500],[pdc1_0; cells_0; mem0; M_0; S_0]);
colors = lines(2);

%% PLOT POPULATION SIZE VS TIME
figure(1); hold on
subplot(4,1,1); hold on;
% yyaxis left; plot(t, y(:, 2), LineWidth=2);
yyaxis left; plot(t, y(:, 2) + y(:, 4), LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
title("Total population size vs. time");

subplot(4,1,2);  hold on;
yyaxis left; plot(t, y(:, 2), LineWidth=2);
title("Number of eng. cells vs time")
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

subplot(4,1,3);  hold on;
yyaxis left; plot(t, y(:, 4), LineWidth=2);
title("Number of cheater cells vs time")
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

subplot(4,1,4);  hold on;
yyaxis left; plot(t, y(:, 1), LineWidth=2); colororder(colors);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
title("Amount of PDC1 vs time");

%% PLOT OVERALL GROWTH RATE VS TIME
total_growth_rate = zeros(length(t), 1);
for i = 1:length(t)
    dxdt = population_odes(t(i), y(i,:)');
    total_growth_rate(i) = dxdt(2) + dxdt(4);
end
figure; plot(t, total_growth_rate, LineWidth=2); title("Population growth rate (cells/hr) vs time")
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

%% PLOT EXP GROWTH TIME CONSTANT VS TIME (ln(N2) - ln(N1))/(t2 - t1)
N = y(:, 2) + y(:, 4);
mu = diff(log(N)) ./ diff(t);
t_mid = (t(1:end-1) + t(2:end)) / 2; % midpoints for plotting

figure;
plot(t_mid, mu, LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
xlabel('Time (h)');
ylabel('\mu (h^{-1})');
title('Population growth rate vs time');

t_mid_2 = (t_mid(1:end-1) + t_mid(2:end)) / 2;
figure; plot(t_mid_2, diff(mu) ./ diff(t_mid), LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

% %% plot ratio of pdc1 to cell count vs. time
% figure;plot(t, y(:,1) ./ y(:, 2), LineWidth=2); title("PDC1:cell count vs. time")