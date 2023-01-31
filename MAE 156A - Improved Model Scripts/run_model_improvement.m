% run_model_improvement.m
%
% Runs model improvement code

%% Clear and set all variables
clc, clear all, close all
% Variable parameters. This function aims to identify improved parameters

% starting guesses to variable parameters
param_var.mu = 0.16;       % coefficient of friction
param_var.cd = 0;           % drag coefficient
param_var.jm = 5e-7;      % motor rotor inertia in [Kg*m^2]

% Fixed parameters
param_fixed.ngear = 4.4;                    % gear ratio
param_fixed.wn = 8200*2*pi/60;              % no load motor speed at 100% PWM in rad/s (converted from rpm)
param_fixed.trq_stall = 0.17*0.0980665;     % motor stall torque at 100% PWM in Nm (converted from kg cm)
param_fixed.nut_mass = 3.20e-3;           % {kg} mass of a single nut
param_fixed.bolt_mass = 7.74e-3;          % {kg} mass of a single bolt
param_fixed.fw_mass = 5.53e-2;              % {kg} flywheel mass
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
config_ar(iconfig).name = ('Config 2 at 75% PWM');
config_ar(iconfig).nut_ar = [1 0 1 0 1 0 1 0];
config_ar(iconfig).pwm = 75;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config22.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 3 at 50% PWM');
config_ar(iconfig).nut_ar = [2 2 2 2 2 2 2 2];
config_ar(iconfig).pwm = 50;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config31.txt'; % Filename
% next configuration
iconfig = iconfig + 1;
config_ar(iconfig).name = ('Config 3 at 100% PWM');
config_ar(iconfig).nut_ar = [2 2 2 2 2 2 2 2];
config_ar(iconfig).pwm = 100;
config_ar(iconfig).exp_data_filename = 'EncoderData_Config33.txt'; % Filename
% set total configurations
nconfig = iconfig;

%% Loop to find optimal coefeciant of friction, param_var.mu
% The assignment will required nested loops to find optimzal values for mu, Jm, and Cd

mu_low_opt = 0.0;
mu_high_opt = 1;
mu_step_opt = 0.25;
mu_vec = mu_low_opt:mu_step_opt:mu_high_opt;
mu_len = length(mu_vec);

Cd_low_opt = 0.0;
Cd_high_opt = .2;
Cd_step_opt = 0.025;
Cd_vec = Cd_low_opt: Cd_high_opt: Cd_step_opt;
Cd_len = length(Cd_vec);

tic
iopt = 1; % mu optimization counter
iopt2 = 1; % Cd optimization counter

for Cd = Cd_low_opt:Cd_step_opt:Cd_high_opt
    param_var.cd = Cd;
    for mu = mu_low_opt:mu_step_opt:mu_high_opt
        disp([Cd mu])
        param_var.mu = mu;
        iconfig = 1;
        for iconfig = 1:nconfig
            config = config_ar(iconfig);
            compare_sim_exp  % run and compare experimental to simulation
        end
        err_metric_mean_opt_ar(iopt2, iopt) = mean(err_metric_ar);
        mu_opt_ar(iopt) = mu;
        
        iopt = iopt + 1;
    end
    Cd_opt_ar(iopt2) = Cd;

    iopt2 = iopt2 + 1;
end
toc

min_err_metric_mean = min(err_metric_mean_opt_ar(err_metric_mean_opt_ar>0));

[r, c] = find(err_metric_mean_opt_ar == min_err_metric_mean, 1, 'last');

fprintf('Most Optimal case corresponds to: \n');
fprintf('Cd of %.4f \n', Cd_vec(r));
fprintf('Mu of %.4f \n \n', mu_vec(c-mu_len*(r-1)));

% figure(4) 
% plot(mu_opt_ar,err_metric_mean_opt_ar)
% xlabel('Coefficient of Friction - \mu');
% ylabel('Mean Error Metric');
% title('Error Metric as a Function of Friction Coefficient');

% output_table. Note use of (:) converts row vectors to column vectors
output_table = [tr_exp_ar(:) wterm_exp_ar(:) tr_sim_ar(:) wterm_sim_ar(:) tr_err_ar(:) wterm_err_ar(:) err_metric_ar(:)];
writematrix(output_table,'Output Table.csv')





