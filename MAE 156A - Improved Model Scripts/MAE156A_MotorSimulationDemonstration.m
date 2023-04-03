%% Optimized ODE45 Motor Model

%MAE 156A Motor Modeling

clear; clc; close all

%% Receive User Inputs

    config.nut_ar = input('Please input nut configuration: '); % Receive nut configuration
    config.pwm = input('Please input PWM value: '); % Receive pwm value
    filename = input('Please input filename: ', 's'); % Receive filename
    param_var.jm = 1.3e-6; % Motor shaft inertia [kg*m^2]
    param_var.mu = 0.18; % Coefficient of Friction
    param_var.cd = 1.3; % Coefficient of Drag
    param_var.tau_f = 0; % Motor friction [Nm]
    
    param_fixed.ngear = 4.4;                    % Gear ratio
    param_fixed.wn = 8150*2*pi/60;              % No load motor speed at 100% PWM in [rad/s]
    param_fixed.trq_stall = 0.0169;             % Motor stall torque at 100% PWM [Nm]
    param_fixed.nut_mass = 3.20e-3;             % Mass of a single nut [kg]
    param_fixed.bolt_mass = 7.74e-3;            % Mass of a single bolt [kg]
    param_fixed.fw_mass = 5.53e-2;              % Flywheel mass [kg]
    param_fixed.hub_mass = 3.84e-3;             % Hub mass [kg]
    param_fixed.r_disk = 0.057;                 % Radius of disk [m]
    param_fixed.r_dis = 0.0508;                 % Distance of nut from center of flywheel [m]

%% Calculate Flywheel Parameters

    [param_var.j_eff, param_var.mfw, param_var.Tau_f] = flywheel_mass_prop(config, param_var, param_fixed);

%% Open Collecteded Data

    collectedData = readmatrix(filename);
    time_exp = collectedData(1:end,1)./1e6;    % [s]
    position_exp = collectedData(1:end,2);    % [counts]
    npt = length(time_exp);

%% Calculate Velocity

    for ipt = 1:npt-1 
        vel(ipt) = (position_exp(ipt+1) - position_exp(ipt)) / ...
            (time_exp(ipt+1) - time_exp(ipt));
    end

%% Calculate Filtered Velocity

    velFilteredRealTime(1) = vel(1);
    alpha = 0.04;  % ranges from 0 to 1
    for ipt = 2:length(vel)
        velFilteredRealTime(ipt) = alpha*vel(ipt) + (1-alpha)*velFilteredRealTime(ipt-1);
    end

%% Plot Experimental Filtered Velocity
    
    figure(1)
    plot(time_exp(2:end),velFilteredRealTime*60/48, 'Color',[0.4940 0.1840 0.5560])
    hold on

%% Motor Simulation and Error Metrics

    [w_sim_arOpt, t_sim_arOpt, tr_simOpt, wterm_simOpt] = motor_sim_ODE45(config, param_var, param_fixed);
    
    [tr_exp, wterm_exp] = find_metrics(velFilteredRealTime*60/48,time_exp(2:end));
    figure(1)
    xline(tr_exp);
    yline(wterm_exp, 'b--');
    
    w_term_err = abs(wterm_exp-wterm_simOpt)/wterm_exp; % percent error in terminal velocity
    tr_err = abs(tr_exp-tr_simOpt)/tr_exp; % percent error in rise time
    err_metric = tr_err + 4*w_term_err; % error metric

    figure(1)
    title_line1 = ['Experimental Velocity (Filtered) vs Optimal jm, mu Optimal cd Simulated Velocity (ODE45) [RPM]'];
    title_line2 = ['Wterm error = ' num2str(w_term_err) '%, tr error = ' num2str(tr_err) '%, Error metric = ' num2str(err_metric)]; % display terminal velocoty with 4 significant figures  
    title_line3 = [' Jm = ' num2str(param_var.jm) ' (kgm^2), Cd = ' num2str(param_var.cd) ', mu = ' num2str(param_var.mu)];
    title({title_line1;title_line2;title_line3})

    xlabel('Time [s]');
    ylabel('Angular Velocity [rpm]');
    legend('Experimental','Optimal jm, mu Optimal cd Simulation','Location','best')