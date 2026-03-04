function dxdt = population_odes(t, x)
    %########## PARAMETERS ##########    
    g_e = .15; d_e = .1;
    g_m = .4; d_m = .1;
    m = @(t) spike(t) + 0.5 * spike(t - 20);
    P_max = 100000; % TODO: vary this with time (P_eff)
    %################################

    e = x(1);
    c = x(2);

    logistic = 1 - (e+c) / P_max;

    dedt = ((g_e - d_e) * e - m(t) * e) * logistic;
    dcdt = ((g_m - d_m) * c + m(t) * e) * logistic;

    dxdt = [dedt; dcdt];
end

function y = spike(t) 
    %########## PARAMETERS ##########
    % r_min; r_max;
    % t0; tf;
    r_min = 5e-7; r_max = 0.9;
    t0 = 10; tf = 15;
    %################################

    slope = (r_max - r_min) / (0.5 * (tf - t0));

    if (t < t0)
        y = r_min;
    elseif (t < (tf + t0)/2)
        y = slope * (t - t0) + r_min;
    elseif (t < tf)
        y = -1 * slope * (t - (tf + t0)/2) + r_max;
    else
        y = r_min;
    end
end

function y = const(t)
    y = 5e-3;
end

%########## INITIAL CONDITIONS ##########
e0 = 1200; c0 = 0;
%########################################

[t,y] = ode45(@population_odes,[0 200],[e0; c0]);

figure;
subplot(3,1,1); plot(t, y, LineWidth=2);
title("Subpopulation size vs. time"); legend("Engineered", "Mutant");

subplot(3,1,2); plot(t, sum(y, 2), LineWidth=2);
title("Total population size vs. time");

subplot(3,1,3); fplot(@(t) spike(t), [0 200], LineWidth=2);
title("Mutation rate vs time");
%%
fplot(@(t) spike(t) + 0.5 * spike(t - 20), [0 200])
