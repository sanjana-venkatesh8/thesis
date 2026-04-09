global light_exp_times light_exp_lengths;

light_exp_times = [23; 59; 98]; %0:30:200;
light_exp_lengths = [15; 11; 24]; %0 * ones(1, length(light_exp_times));

global k_delay;
global r_max r_min;
global P_max;
global g_pdc1 r_cons R_p;
global g_s d_s R_s;
global g_N d_N;
global g_c d_C;
global k_leak;
global q_pdc q_N q_M;

k_delay = .01; %[done]
r_max = 5e-3; r_min = 1e-6; %[done]
P_max = 8e5; %[done]

g_pdc1 = .5; r_cons = .2; R_p = 80; %[done]
g_s = .05; d_s = 0.01; R_s = 3; % do with Sp Batch
g_N = 0.5; d_N = 0.01; % do with ON/OFF
d_C = 0.05; %[done]
g_c = 0.15; %[done]
k_leak = 1e-9; %[done]
q_pdc = 1e-5; q_N = 1e-5; q_M = 1000; 
% q_pdc = 0; q_N = 1e-5; q_M = 1000; % do with Sp Batch

function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########  
    global k_delay;
    global r_max r_min;
    global P_max;
    global g_pdc1 r_cons R_p;
    global g_s d_s R_s;
    global g_N d_N;
    global g_c d_C;
    global k_leak;
    global q_pdc q_N q_M; % how much does adding SpSwi6/Clr4 attenuate cell processes
    %################################

    pdc1 = x(1);
    s = x(2);
    N = x(3);
    C = x(4);
    mem = x(5);
    rm = x(6);

    cell_logistic = 1;%(1 - (N + C)/P_max);

    dNdt = ((g_N) / (g_N*q_N*s + 1) * pdc1 - d_N * N) * cell_logistic;
    dsdt = (g_s - d_s) * N * (1 - s/(N * R_s));
    d_pdc1_dt = ((light(t)*g_pdc1) / (light(t)*g_pdc1*q_pdc*s + 1) * N) * (1 - pdc1/(N * R_p)) - r_cons * dNdt;
    dCdt = ((rm) / (rm * q_M * s + 1) * N - d_C * C + g_c * C) * cell_logistic;
    d_mem_dt = (1-light(t)) - k_leak*mem;
    d_rm_dt = (k_delay*mem)*(r_max*(1-light(t)) - rm + r_min); % cannot use this to show slowed mutation rate
    d_rmEff_dt = (rm) / (rm * q_M * s + 1);

    dxdt = [d_pdc1_dt; dsdt; dNdt; dCdt; d_mem_dt; d_rm_dt; d_rmEff_dt];
    % dxdt = [d_pdc1_dt; dsdt; dNdt; dCdt; d_mem_dt; d_rm_dt; d_rmEff_dt];
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

function plotbg(args)
    arguments
        args.trial {mustBeMember(args.trial,["ctrl","sp"])}
        args.avg logical
        args.batch {mustBeMember(args.batch, ["tiara", "sergio"])}
    end

    trial = args.trial;
    avg = args.avg;
    batch = args.batch;

    if batch == "sergio"
        if trial == "ctrl"
            data = readmatrix("data_vis/ALL_dark.xlsx", Sheet="139-control");
        else
            data = readmatrix("data_vis/ALL_dark.xlsx", Sheet="162(SpSwi6-SpClr4)");
        end
    
        gcf; hold on;
    
        if avg
            plot(data(:,1), mean(data(:,2:end), 2, 'omitnan'), LineWidth = 2, Color=[0.5 0 1 0.1],Marker="none");
        else
            plot(data(:,1), data(:,2:end), LineWidth = 2, Color=[0.5 0 1 0.1],Marker="none");
        end
    else
        if trial == "ctrl"
            times = readmatrix("data_vis/cc0_times.csv");
            data = readmatrix("data_vis/cc0_data.csv");
        else
            times = readmatrix("data_vis/cc1_times.csv");
            data = readmatrix("data_vis/cc1_data.csv");
        end

        gcf; hold on;

        if avg
            plot(mean(times,2, 'omitnan'), mean(data, 2, 'omitnan'), LineWidth = 2, Color=[0 0.5 1 0.1],Marker="none");
        else
            hold on;
            for i = 1:size(times, 2)
                plot(times(:,i), data(:,i), LineWidth = 2, Color=[0.3 0.6 0.45 0.05],Marker="none");
            end
        end
    end
    hold off;
end
%%
%########## INITIAL CONDITIONS ##########
pdc1_0 = 35000;%35000; 
s_0 = 0;
N_0 = 0.01e6; C_0 = 0;
mem0 = 0; rm0 = 1e-3;
%########################################

[t,y] = ode45(@population_odes,[0 480],[pdc1_0; s_0; N_0; C_0; mem0; rm0; rm0]);

titles = ["PDC concentration vs time", "Chromatin binding protein concentration vs. time", "N cell population vs time", "Cheater population vs. time", "Total time in darkness vs. time", "Mutation rate vs. time", "Effective mutation rate vs. time"];
figure(1); hold on;
for i = 1:size(y, 2)
    subplot(size(y, 2), 1, i); hold on;
    plot(t, y(:, i), LineWidth=2); title(titles(i));
    xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
end

figure(2);
yyaxis right; plotbg(trial="sp", avg=false, batch="tiara");
plotbg(trial="sp", avg=false, batch="sergio");
hold on; yyaxis left;
plot(t, y(:, 3) + y(:, 4), LineWidth=2, Marker='none'); title("total cell population vs. time")
plot(t, y(:, 3), LineWidth=2, LineStyle="--", DisplayName="Engineered cells");
plot(t, y(:, 4), LineWidth=2, LineStyle='--', DisplayName="Cheaters");
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);

ax = gca;
ax.YAxis(1).Limits = [N_0 P_max]; % Set identical limits
ax.YAxis(2).Limits = [0.01 7];

%% PLOT EXP GROWTH TIME CONSTANT VS TIME (ln(N2) - ln(N1))/(t2 - t1)
N = y(:, 3);
C = y(:, 4);
total = N + C;
mu_N = diff(log(N)) ./ diff(t);
mu_C = diff(log(C)) ./ diff(t);
mu_total = diff(log(total)) ./ diff(t);
t_mid = (t(1:end-1) + t(2:end)) / 2; % midpoints for plotting

figure(3); hold on;
plot(t_mid, mu_total, LineWidth=2);
plot(t_mid, mu_C, LineWidth=2, LineStyle='--', DisplayName="C");
plot(t_mid, mu_N, LineWidth=2, LineStyle='--', DisplayName="N");
ylim([0 0.5])
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
xlabel('Time (h)');
ylabel('\mu (h^{-1})');
title('Population growth rate vs time');

% %% PLOT OVERALL GROWTH RATE VS TIME
% total_growth_rate = zeros(length(t), 1);
% for i = 1:length(t)
%     dxdt = population_odes(t(i), y(i,:)');
%     total_growth_rate(i) = dxdt(3) + dxdt(4);
% end
% % figure; 
% hold on; plot(t, total_growth_rate, LineWidth=2); title("Population growth rate (cells/hr) vs time")
% xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
