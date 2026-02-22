%% DESCRIPTION
% Mutations are binned into k types, but each mutant cell has its own
% fitness s - the rationale is that there are many unique mutations, but
% they affect the cell in a limited number of ways.

% Another piece of rationale is from Eigen et al. - selection pressures act
% on groups of mutations in similar ways.

%% MODEL PARAMETERS
% initial population size (10mL of cells at OD0.01 (120 cells/mL)
E0 = 1200;

% Carrying capacity (for batch culture of only E cells)
K    = 100;      % tune as needed

% Base division rates
g_E  = 0.3;     % (h^-1) engineered (burdened)
g_C  = 1.0;     % (h^-1) cheater (relieved burden, fitter)

% Death rate
r_death  = 0.03;    % (h^-1)

% Mutation *probability* - for simplicity, this is the same for E and C cells
% TODO: time-evolving mutation rate (for control strain only)
r_mutation = 0.04;

% DFE parameters
% Engineered: many beneficial mutations (loss of function / lower burden)
f_b = 0.5;   f_n = 0.2;  f_d = 0.3;
lambda_b = 1.5;   % larger mean beneficial effect (1/lambda)
lambda_d = 8;

% TIME/REP PARAMETERS
Tmax     = 200;  % max time (hours)
num_reps = 100;

%% PLOT DFE
% TODO

%% GILLESPIE SIMULATION

n_mutation_classes = 5;

all_traj = cell(num_reps, 1); 

figure; hold on; % DEBUG

for rep = 1:num_reps
    tic
    
    t  = 0;
    counts = zeros(n_mutation_classes, 1);  % tracks # of cells in each mutation class
    counts(1) = E0;                         % the first mutation class is the engineered class
    mean_s  = zeros(n_mutation_classes, 1); % average fitness per class

    % TODO: assign some random s
    mean_s(1) = 0.3;

    T_traj = t;
    N_traj = E0; % initial population size = all engineered cells

    while t < Tmax && sum(counts) > 0
        N = sum(counts);

        alpha  = 0.1;
        K_eff  = K * (1 + alpha*mean_s); % set K_eff based on avg fitness
        logistic = max(0, 1 - N / K_eff);

        % Average division, death, mutation propensities for each mutation class
        a_div = mean_s .* logistic.'; % dimensions wrong
        a_death = ones(n_mutation_classes, 1) .* r_death;
        a_mut = ones(n_mutation_classes, 1) .* r_mutation;
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

        % Pick mutation class
        u2 = rand * a0;
        cumul = cumsum(a_i);
        idx = find(u2 < cumul, 1, 'first');

        % Pick which event (duplication, death, mutation)
        i = rand;
        if i < a_div(idx) / a_i(idx)
            % Division event:
            % 1. offspring is of same mutation class as parent
            % 2. parent class fitness does not change
            counts(idx) = counts(idx) + 1;

        elseif i < (a_div(idx) + a_death(idx)) / a_i(idx)
            % Death event
            counts(idx) = max(0, counts(idx) - 1);

        else
            % Mutation event
            % 1. offspring is of a new mutation class
            % 2. New mutation class fitness is changed by an amount dS

            mutation_class = randi([2 n_mutation_classes]);

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

            % Add offspring and update fitness
            mean_s(mutation_class) = (mean_s(mutation_class) * counts(mutation_class) + dS) / (counts(mutation_class) + 1);
            counts(mutation_class) = counts(mutation_class) + 1;
        end

        % Record trajectory
        T_traj(end+1, 1) = t; %#ok<AGROW>
        N_traj(end+1, 1) = sum(counts);
    end

    all_traj{rep} = struct('t', T_traj, 'N', N_traj);
    
    fprintf("Rep %d completed in ", rep);
    toc

end

%% Plot total population trajectories
figure(1); hold on;
% figure(3); hold on; % DEBUG: plot avg fitness v time
% figure(4); hold on; % DEBUG: plot avg mutant fitness v time
% figure(4); hold on; % DEBUG: plot avg eng fitness v time
cols = lines(num_reps);

for rep = 1:num_reps
    tr = all_traj{rep};
    figure(1); plot(tr.t, tr.N, 'Color', cols(rep,:), 'LineWidth', 1.3);
end
figure(1); 
xlabel('Time (hours)');
ylabel('Total population size');
title('Cell count vs. time');
grid on;

% figure(3); 
% xlabel('Time (hours)');
% ylabel('Average population fitness');
% title('Average population fitness vs. time');
% grid on;
% 
% figure(4); title("Average mutant fitness")
% figure(5); title("Average eng fitness")
