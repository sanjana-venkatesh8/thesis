global light_exp_times light_exp_lengths;
light_exp_times = 0:30:200;%[10, 30, 50, 70, 90];
light_exp_lengths = 3 * ones(1, length(light_exp_times));

function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########    
    g_cells = 0.23; r_death = 0.15; cells_max = 100000;
    g_pdc1 = 3;
    r_pdcCons = 2;
    R = 8;
    %################################

    pdc1 = x(1);
    cells = x(2);

    dcellsdt = (g_cells*pdc1-r_death*cells);%*(1-cells/cells_max);
    dpdc1dt = (light(t)*g_pdc1*cells - dcellsdt * r_pdcCons);% * (1 - pdc1/(cells*R));

    dxdt = [dpdc1dt; dcellsdt];
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
pdc1_0 = 0.15; cells_0 = 1200;
%########################################

[t,y] = ode45(@population_odes,[0 200],[pdc1_0; cells_0]);
colors = lines(2);

%% PLOT POPULATION SIZE VS TIME
figure(1); hold on
subplot(2,1,1); hold on;
yyaxis left; plot(t, y(:, 2), LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
title("Total population size vs. time");

subplot(2,1,2);  hold on;
yyaxis left; plot(t, y(:, 1), LineWidth=2); colororder(colors);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
title("Amount of PDC1 vs time");

% %% PLOT OVERALL GROWTH RATE VS TIME
% total_growth_rate = zeros(length(t), 1);
% for i = 1:length(t)
%     dxdt = population_odes(t(i), y(i,:)');
%     total_growth_rate(i) = dxdt(2);
% end
% figure; plot(t, total_growth_rate, LineWidth=2)

%% PLOT EXP GROWTH TIME CONSTANT VS TIME (ln(N2) - ln(N1))/(t2 - t1)
N = y(:, 2);
mu = diff(log(N)) ./ diff(t);
t_mid = (t(1:end-1) + t(2:end)) / 2; % midpoints for plotting

figure;
plot(t_mid, mu, LineWidth=2);
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
xlabel('Time (h)');
ylabel('\mu (h^{-1})');
title('Population growth rate vs time');

%% plot ratio of pdc1 to cell count vs. time
figure;plot(t, y(:,1) ./ y(:, 2), LineWidth=2)