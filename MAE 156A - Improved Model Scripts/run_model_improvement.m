% run_model_improvement.m
%
% Runs model improvement code

%% Clear and set all variables
clc, clear all, close all
% Variable parameters. This function aims to identify improved parameters

% starting guesses to variable parameters
param_var.mu = 0.18;        % coefficient of friction
param_var.cd = 1.3;           % drag coefficient
param_var.jm = 1.3e-6;        % motor rotor inertia in [Kg*m^2]
param_var.Tau_f = 0;        % motor friction [N*m]

% Fixed parameters
param_fixed.ngear = 4.4;                    % gear ratio
param_fixed.wn = 8150*2*pi/60;              % no load motor speed at 100% PWM in rad/s (converted from rpm)
param_fixed.trq_stall = 0.0169;             % motor stall torque at 100% PWM in Nm
param_fixed.nut_mass = 3.20e-3;             % {kg} mass of a single nut
param_fixed.bolt_mass = 7.74e-3;            % {kg} mass of a single bolt
param_fixed.fw_mass = 4.4e-2;              % {kg} flywheel mass
param_fixed.hub_mass = 3.84e-3;             % {kg} hub mass
param_fixed.r_disk = 0.057;                 % radius of disk [m]
param_fixed.r_dis = 0.0508;                 % distance of nut from center[m]
param_fixed.gearbox_fric_factor = 6;        % factor that multiplies coef of friction

% Configuration variables
iconfig=1;          % Index number
config_ar(iconfig).name = ('Config 1 at 50% PWM');
config_ar(iconfig).nut_ar = [1 0 0 0 1 0 0 0];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config11.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 1 at 75% PWM');
config_ar(iconfig).nut_ar = [1 0 0 0 1 0 0 0];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config12.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 1 at 100% PWM');
config_ar(iconfig).nut_ar = [1 0 0 0 1 0 0 0];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config13.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 2 at 50% PWM');
config_ar(iconfig).nut_ar = [1 0 1 0 1 0 1 0];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config21.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 2 at 75% PWM');
config_ar(iconfig).nut_ar = [1 0 1 0 1 0 1 0];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config22.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 2 at 100% PWM');
config_ar(iconfig).nut_ar = [1 0 1 0 1 0 1 0];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config23.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 3 at 50% PWM');
config_ar(iconfig).nut_ar = [2 2 2 2 2 2 2 2];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config31.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 3 at 75% PWM');
config_ar(iconfig).nut_ar = [2 2 2 2 2 2 2 2];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config32.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 3 at 100% PWM');
config_ar(iconfig).nut_ar = [2 2 2 2 2 2 2 2];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config33.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 4 at 50% PWM');
config_ar(iconfig).nut_ar = [3 0 0 0 3 0 0 0];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config41.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 4 at 75% PWM');
config_ar(iconfig).nut_ar = [3 0 0 0 3 0 0 0];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config42.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 4 at 100% PWM');
config_ar(iconfig).nut_ar = [3 0 0 0 3 0 0 0];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config43.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 5 at 50% PWM');
config_ar(iconfig).nut_ar = [3 0 3 0 3 0 3 0];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config51.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 5 at 75% PWM');
config_ar(iconfig).nut_ar = [3 0 3 0 3 0 3 0];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config52.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 5 at 100% PWM');
config_ar(iconfig).nut_ar = [3 0 3 0 3 0 3 0];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config53.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 6 at 50% PWM');
config_ar(iconfig).nut_ar = [3 3 3 3 3 3 3 3];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config61.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 6 at 75% PWM');
config_ar(iconfig).nut_ar = [3 3 3 3 3 3 3 3];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config62.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 6 at 100% PWM');
config_ar(iconfig).nut_ar = [3 3 3 3 3 3 3 3];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config63.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 7 at 50% PWM');
config_ar(iconfig).nut_ar = [1 1 1 1 1 1 1 1];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config71.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 7 at 75% PWM');
config_ar(iconfig).nut_ar = [1 1 1 1 1 1 1 1];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config72.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 7 at 100% PWM');
config_ar(iconfig).nut_ar = [1 1 1 1 1 1 1 1];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config73.txt'; % Filename
% set total configurations
nconfig = iconfig;

%% Loop to find optimal coefeciant of friction, param_var.mu
% The assignment will required nested loops to find optimzal values for mu, Jm, and Cd

Opt_case = 'mu'; % Variable that is being optimized

if contains(Opt_case, 'mu')

    mu_low_opt = 0.17;
    mu_high_opt = 0.19;
    mu_step_opt = 0.001;
    mu_vec = mu_low_opt:mu_step_opt:mu_high_opt;
    mu_len = length(mu_vec);
    
    tic
    iopt = 1; % mu optimization counter

    for mu = mu_low_opt:mu_step_opt:mu_high_opt
    disp(mu)
    param_var.mu = mu;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    mu_opt_ar(iopt) = mu;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(mu_opt_ar,err_metric_mean_opt_ar)
    xlabel('Coefficient of Friction - \mu');
    ylabel('Mean Error Metric');
    title('Error Metric as a Function of Friction Coefficient');

elseif contains(Opt_case, 'cd')

    cd_low_opt = 1.0;
    cd_high_opt = 1.4;
    cd_step_opt = 0.05;
    cd_vec = cd_low_opt: cd_high_opt: cd_step_opt;
    cd_len = length(cd_vec);

    tic
    iopt = 1; % cd optimization counter

    for cd = cd_low_opt:cd_step_opt:cd_high_opt
    disp(cd)
    param_var.cd = cd;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    cd_opt_ar(iopt) = cd;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(cd_opt_ar,err_metric_mean_opt_ar)
    xlabel('Coefficient of Drag - cd');
    ylabel('Mean Error Metric');
    title('Error Metric as a Function of Drag Coefficient');

elseif contains(Opt_case, 'jm')

    jm_low_opt = 1.2e-6;
    jm_high_opt = 1.4e-6;
    jm_step_opt = 5e-8;
    jm_vec = jm_low_opt: jm_high_opt: jm_step_opt;
    jm_len = length(jm_vec);

    tic
    iopt = 1; % cd optimization counter

    for jm = jm_low_opt:jm_step_opt:jm_high_opt
    disp(jm)
    param_var.jm = jm;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    jm_opt_ar(iopt) = jm;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(jm_opt_ar,err_metric_mean_opt_ar)
    xlabel('Coefficient of Shaft Inertia - jm');
    ylabel('Mean Error Metric');
    title('Error Metric as a Function of Shaft Inertia');

elseif contains(Opt_case, 'wn')

    tr_stall_low_opt = 8100;
    trq_stall_high_opt = 8200;
    trq_stall_step_opt = 10;
    trq_stall_vec = tr_stall_low_opt: trq_stall_high_opt: trq_stall_step_opt;
    trq_stall_len = length(trq_stall_vec);

    tic
    iopt = 1; % cd optimization counter

    for trq_stall = tr_stall_low_opt:trq_stall_step_opt:trq_stall_high_opt
    disp(trq_stall)
    param_fixed.wn = trq_stall*2*pi/60;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    trq_stall_opt_ar(iopt) = trq_stall;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(trq_stall_opt_ar,err_metric_mean_opt_ar)
    xlabel('No Load speed [RPM]');
    ylabel('Mean Error Metric');
    title('Error Metric as a No Load Speed');

elseif contains(Opt_case, 'trq_stall')

    tr_stall_low_opt = 0.016;
    trq_stall_high_opt = 0.018;
    trq_stall_step_opt = 0.0001;
    trq_stall_vec = tr_stall_low_opt: trq_stall_high_opt: trq_stall_step_opt;
    trq_stall_len = length(trq_stall_vec);

    tic
    iopt = 1; % cd optimization counter

    for trq_stall = tr_stall_low_opt:trq_stall_step_opt:trq_stall_high_opt
    disp(trq_stall)
    param_fixed.trq_stall = trq_stall;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    trq_stall_opt_ar(iopt) = trq_stall;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(trq_stall_opt_ar,err_metric_mean_opt_ar)
    xlabel('Stall Torque [Nm]');
    ylabel('Mean Error Metric');
    title('Error Metric as a Function of Stall Torque');

elseif contains(Opt_case, 'fw_mass')

    fw_mass_low_opt = 3.5e-2;
    fw_mass_high_opt = 4.5e-2;
    fw_mass_step_opt = 5e-4;
    fw_mass_vec = fw_mass_low_opt: fw_mass_high_opt: fw_mass_step_opt;
    fw_mass_len = length(fw_mass_vec);

    tic
    iopt = 1; % cd optimization counter

    for fw_mass = fw_mass_low_opt:fw_mass_step_opt:fw_mass_high_opt
    disp(fw_mass)
    param_fixed.fw_mass = fw_mass;
    iconfig = 1;
        for iconfig = 1:nconfig
            fclose all;
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
    err_metric_mean_opt_ar(iopt) = mean(err_metric_ar);
    fw_mass_opt_ar(iopt) = fw_mass;
    
    iopt = iopt + 1;
    end
    toc

    figure(4) 
    plot(fw_mass_opt_ar,err_metric_mean_opt_ar)
    xlabel('Flywheel Mass [kg]');
    ylabel('Mean Error Metric');
    title('Error Metric as a Function of Flywheel Mass');

elseif contains(Opt_case, 'load_var')

end






% min_err_metric_mean = min(err_metric_mean_opt_ar(err_metric_mean_opt_ar>0));
% 
% [r, c] = find(err_metric_mean_opt_ar == min_err_metric_mean, 1, 'last');

% fprintf('Most Optimal case corresponds to: \n');
% fprintf('Cd of %.4f \n', Cd_vec(r));
% fprintf('Mu of %.4f \n \n', mu_vec(c-mu_len*(r-1)));


% output_table. Note use of (:) converts row vectors to column vectors
output_table = [tr_exp_ar(:) wterm_exp_ar(:) tr_sim_ar(:) wterm_sim_ar(:) tr_err_ar(:) wterm_err_ar(:) err_metric_ar(:)];
writematrix(output_table,'Output Table.csv')





