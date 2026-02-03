%% MODEL PARAMETERS
% initial population size (10mL of cells at OD0.01 (120 cells/mL)
E0 = 1200;

% Carrying capacity (for batch culture of only E cells)
K    = E0;      % tune as needed

% Base division rates
g_E  = 0.3;      % engineered (burdened)
g_C  = 1.0;      % cheater (relieved burden, fitter)

% Death rates
d_E  = 0.02;
d_C  = 0.02;

% Mutation probability - for simplicity, this is the same for E and C cells
% TODO: time-evolving mutation rate (for control strain only)
mu = 5e-5;

% DFE parameters
% Engineered: many beneficial mutations (loss of function / lower burden)
f_b = 0.5;   f_n = 0.2;  f_d = 0.3;
lambda_b = 1.5;   % larger mean beneficial effect (1/lambda)
lambda_d = 8;

% Cheater: closer to optimum, beneficials rarer/weaker
f_b_C = 0.01;  f_n_C = 0.39; f_d_C = 0.60;
lambda_b_C = 10;
lambda_d_C = 10;

Tmax     = 200;  % max time (hours)
num_reps = 20;

%% GILLESPIE SIMULATION
% TODO: track the number of E- and C-cells
% TODO: for mutating cheater cells, s sample gets added onto old s

all_traj = cell(num_reps, 1); % DEBUG

figure; hold on;

for rep = 1:num_reps
    tic
    
    t  = 0;
    B  = repmat('E', E0, 1);  % tracks individual cells; 'E' = engineered, 'C' = cheater
    s  = zeros(E0, 1);        % fitness modifiers for each cell (relative to g_E or g_C)

    T_traj = t;
    N_traj = E0; % initial population size = all engineered cells

    while t < Tmax && numel(B) > 0
        isE = (B == 'E');
        isC = ~isE;

        N = numel(B);
        n_E = numel(B(isE));
        n_C = numel(B(isC));

        alpha  = 0.7;
        mean_s = mean(s); % avg fitness of population
        K_eff  = K * (1 + alpha * mean_s); % set K_eff based on avg fitness
        logistic = max(0, 1 - N / K_eff);

        % Division, death, mutation propensities for each cell
        a_div = zeros(N, 1);
        % a_div(isE) = g_E * n_E * (1 + s(isE)) * logistic;
        % a_div(isC) = g_C * n_C * (1 + s(isC)) * logistic;
        a_div(isE) = g_E * (1 + s(isE)) * logistic;
        a_div(isC) = g_C * (1 + s(isC)) * logistic;


        a_death = zeros(N, 1);
        % a_death(isE) = d_E * n_E;
        % a_death(isC) = d_C * n_C;
        a_death(isE) = d_E;
        a_death(isC) = d_C;

        % TODO: LOOK AT THIS CAREFULLY
        a_mut = zeros(N, 1);
        % a_mut(isC) = mu * n_E; % both of these get added to C
        % a_mut(isC) = a_mut(isC) + mu * n_C;
        a_mut(isC) = mu; % both of these get added to C
        a_mut(isC) = a_mut(isC) + mu;

        a_i = a_div + a_death + a_mut;
        a0 = sum(a_i);

        if (a0 <= 0)
            break;
        end

        % Sample time to next event
        u1 = rand;
        tau = -log(u1) / a0;
        t   = t + tau;
        if t > Tmax
            break;
        end

        % Pick cell
        u2 = rand * a0;
        cumul = cumsum(a_i);
        idx = find(u2 < cumul, 1, 'first');

        % Pick which event (duplication, death, mutation)
        i = rand;
        parentB = B(idx);
        parentS = s(idx);
        if i < a_div(idx) / a_i(idx)
            % Division event - offspring has same status (E/C) and fitness
            % as parent
            offB = parentB;
            offS = parentS;

            % Add offspring
            B(end+1,1) = offB;
            s(end+1,1) = offS;
        elseif i > a_div(idx) / a_i(idx) && i < (a_div(idx) + a_death(idx)) / a_i(idx)
            % Death event
            B(idx,:) = [];
            s(idx,:) = [];
        else
            % Mutation event
            offB = 'C';
            % draw from DFE to get dS
            u = rand;
            if u < f_b
                dS = exprnd(1 / lambda_b);
            elseif u < f_b + f_n
                dS = 0;
            else
                % Deleterious: Gamma(k, theta) with mean k*theta ~ 1/lambda_d_E
                k_del   = 2;                        % shape (tunable)
                theta_d = 1 / (lambda_d * k_del); % scale so mean = 1/lambda_d_E
                dS      = -gamrnd(k_del, theta_d);  % negative deleterious effect
            end
            offS = parentS + dS;

            % Add offspring
            B(end+1,1) = offB;
            s(end+1,1) = offS;
        end

        % Record trajectory
        T_traj(end+1, 1) = t; %#ok<AGROW>
        N_traj(end+1, 1) = numel(B);


    end

    all_traj{rep} = struct('t', T_traj, 'N', N_traj);
    
    fprintf("Rep %d completed in ", rep);
    toc

    figure(1); scatter(rep, numel(B(isC))) % DEBUG
end

%% Plot total population trajectories
figure; hold on;
cols = lines(num_reps);
for rep = 1:num_reps
    tr = all_traj{rep};
    plot(tr.t, tr.N, 'Color', cols(rep,:), 'LineWidth', 1.3);
end
xlabel('Time (hours)');
ylabel('Total population size');
title('SSA trajectories: engineered-only start, cheaters emerge over time');
grid on;