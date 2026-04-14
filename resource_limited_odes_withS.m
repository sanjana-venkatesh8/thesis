%% DECLARE PARAMETERS
global light_exp_times light_exp_lengths;

global k_delay; % IMPORTANT
global r_max r_min;
global P_max;
global g_pdc1 d_pdc1 R_p;
global g_s d_s R_s;
global g_N d_N;
global g_c d_C;
global k_leak; % IMPORTANT
global q_pdc q_N q_M; % q_M IMPORTANT

global culture_type;
global strain_type;

%% SET PARAMETER VALUES AND INITIAL CONDITIONS

light_exp_times = 0:20:200; %[0; 38; 70];
light_exp_lengths = 10 * ones(1, length(light_exp_times)); %[23; 21; 28];

k_delay = 100;
r_max = 5e-3; r_min = 1e-6;
P_max = 8e5;

g_pdc1 = .5; d_pdc1 = .1; R_p = 8;
g_s = .05; d_s = 0.01; R_s = 3;
g_N = 0.25; d_N = 0.005;
d_C = 0.005;
g_c = 0.1;
k_leak = 1e-9;
q_pdc = 3e-3; q_N = 5e-5; q_M = 1000;

culture_type = "turbidostat"; % "batch" or "turbidostat"
strain_type = "ctrl"; %"ctrl" or "sp"

%########## INITIAL CONDITIONS ##########
pdc1_0 = 0;%30000;
s_0 = 0;
N_0 = 0.01e6; 
C_0 = 0;
mem0 = 0; rm0 = 1e-3;
%########################################

%% ODE MODEL AND HELPER FUNCTIONS
function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########  
    global k_delay;
    global r_max r_min;
    global P_max;
    global g_pdc1 d_pdc1 R_p;
    global g_s d_s R_s;
    global g_N d_N;
    global g_c d_C;
    global k_leak;
    global q_pdc q_N q_M; % how much does adding SpSwi6/Clr4 attenuate cell processes

    global culture_type; global strain_type;
    %################################

    pdc1 = x(1);
    s = x(2);
    N = x(3);
    C = x(4);
    mem = x(5);
    rm = x(6);

    if culture_type == "batch"
        cell_logistic = (1 - (N + C)/P_max);
    else % if culture_type == 'turbidostat'
        cell_logistic = 1;
    end

    dNdt = ((g_N) / (g_N*q_N*s + 1) * pdc1 - d_N * N) * cell_logistic;
    dsdt = (g_s*N) * (1 - s/(N * R_s)) - d_s*s;
    d_pdc1_dt = ((light(t)*g_pdc1) / (light(t)*g_pdc1*q_pdc*s + 1) * N) * (1 - pdc1/(N * R_p)) - d_pdc1 * pdc1;
    dCdt = ((rm) / (rm * q_M * s + 1) * N - d_C * C + g_c * C) * cell_logistic;
    d_mem_dt = (1-light(t)) - k_leak*mem;
    d_rm_dt = (1/k_delay*mem)*(r_max*(1-light(t)) - rm + r_min); % cannot use this to show slowed mutation rate
    d_rmEff_dt = (rm) / (rm * q_M * s + 1);

    if strain_type == "ctrl"
        % CONTROL STRAIN
        dxdt = [d_pdc1_dt; 0; dNdt; 0; 0; 0; 0];
    else % if strain_type == "sp"
        % CHROMATIN-COMPACTED STRAIN
        dxdt = [d_pdc1_dt; dsdt; dNdt; dCdt; d_mem_dt; d_rm_dt; d_rmEff_dt];
    end
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

function non_exposed_time = time_in_dark(t)
    % Inputs:
    %   t: cheating time (>0)
    % Output:
    %   non_exposed_time: total time NOT exposed in [0, t]
    
    global light_exp_times light_exp_lengths; 

    n = length(light_exp_times);
    exp_starts = light_exp_times(:);  % Ensure column vector
    exp_ends = exp_starts + light_exp_lengths(:);
    
    % Filter and clip exposures overlapping [0,t]
    valid_mask = (exp_starts < t) & (exp_ends > 0);
    if ~any(valid_mask)
        total_exposed = 0;
    else
        valid_starts = max(exp_starts(valid_mask), 0);
        valid_ends = min(exp_ends(valid_mask), t);
        total_exposed = sum(valid_ends - valid_starts);
    end
    
    total_time = t;
    non_exposed_time = total_time - total_exposed;
end

%% LIGHT CYCLE TIME VS DURATION SWEEP
% time_vals = 10:10:100;
% duration_vals = 0:5:100;
% 
% no_light_duration = nan(length(time_vals), length(duration_vals));
% 
% for i = 1:length(time_vals)
%     fprintf("%d: ", i)
%     light_exp_times = 0:time_vals(i):500; 
%     max_duration_ind = length(duration_vals(1:find(duration_vals>time_vals(i), 1, 'first'))) - 1;
%     if max_duration_ind == -1
%         max_duration_ind = length(duration_vals);
%     end
%     for j = 1:max_duration_ind
%         light_exp_lengths = duration_vals(j) * ones(1, length(light_exp_times)); 
%         [t,y] = ode45(@population_odes,[0 500],[pdc1_0; s_0; N_0; C_0; mem0; rm0; rm0]);
%         N = y(:,3);
%         C = y(:,4);
% 
%         t_esc = t(find(C > N, 1, 'first'));
%         if isscalar(t_esc)
%             % no_light_duration(i, j) = time_in_dark(t_esc);
%             no_light_duration(i, j) = t_esc;
%         else
%             no_light_duration(i, j) = inf;
%         end
%         fprintf("%d ", j)
%     end
%     fprintf("\n")
% end
% 
% figure; 
% heatmap(duration_vals, time_vals, no_light_duration);
% ylabel("time"); xlabel("duration");
%% SOLVE AND PLOT SOLUTIONS

if culture_type == "batch"
    light_exp_times = 0;
    light_exp_lengths = 0;
end

[t,y] = ode45(@population_odes,[0 200],[pdc1_0; s_0; N_0; C_0; mem0; rm0; rm0]);

titles = ["PDC concentration vs time", "Chromatin binding protein concentration vs. time", "N cell population vs time", "Cheater population vs. time", "Total time in darkness vs. time", "Mutation rate vs. time", "Effective mutation rate vs. time"];
figure(1); hold on;
for i = 1:size(y, 2)
    subplot(size(y, 2), 1, i); hold on;
    plot(t, y(:, i), LineWidth=2); title(titles(i));
    xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
end

% PLOT TOTAL CELL POPULATION VS TIME
figure(2);
yyaxis right; plotbg(trial="ctrl", avg=false, batch="tiara");
plotbg(trial="ctrl", avg=false, batch="sergio");
hold on; yyaxis left;
plot(t, y(:, 3) + y(:, 4), LineWidth=2, Marker='none'); title("total cell population vs. time")
plot(t, y(:, 3), LineWidth=2, LineStyle="--", DisplayName="Engineered cells");
plot(t, y(:, 4), LineWidth=2, LineStyle='--', DisplayName="Cheaters");
xregion(light_exp_times, light_exp_times+light_exp_lengths, FaceColor=[100/255, 225/255, 240/255], FaceAlpha=0.2);
yyaxis left; 
ylabel("Number of cells (Simulation)"); 
yyaxis right; ylabel("OD_{600} (experiment)")
xlabel("Time (hours)");

ylim([N_0 P_max])
ax = gca;
ax.YAxis(1).Limits = [N_0 P_max]; % Set identical limits for simulation and data
ax.YAxis(2).Limits = [0.01 7];

% PLOT EXP GROWTH TIME CONSTANT VS TIME (ln(N2) - ln(N1))/(t2 - t1)
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