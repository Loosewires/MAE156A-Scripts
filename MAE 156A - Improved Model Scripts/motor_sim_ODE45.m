function [w_sim_ar, t_sim_ar, tr_sim, wterm_sim] = motor_sim_ODE45(config, param_var, param_fixed)
% ODE45 method simulation of motor

pl_flag = false; % true (1) or false (0) to use to turn plotting on or off

% Motor Parameters
w_nl = param_fixed.wn*config.pwm/100; % No load speed in [rad/s] at PWM percentatge
tau_s =  param_fixed.trq_stall*config.pwm/100; % Stall torque in [N-m] at PWM percentage
I_r = param_var.j_eff;  % effective rotational inertia reflected to motor in [kg*m^2]
rs = .002; % Radius of output shaft [m]
g = 9.8; % Acceleration due to gravity [m/s^2]
Lc = .0099; % Length from load to left part [m]
Lb = .0066; % length from left part to right part [m]
rho_air = 1.204; % Density of air at 20C [kg/m^3]
A_cross = .005*.01; % Cross sectional area of bolt/nuts [m^2]
r_bolt = .002; % Radius of bolt [m]

% Simulation Parameters are hard coded for convinence
tspan = [0 10]; % Time for simulation

% Initial conditions
w0 = 0; % Initial angular velocity [rad/s]

% Calculate Friction
Tau_f = param_var.mu*rs*g*param_var.mfw/param_fixed.ngear * (2*(Lc+Lb)/Lb - 1);

% Calculate Constant part of Drag
Tau_d = 1/2*rho_air*A_cross*r_bolt^3/param_fixed.ngear^2;

% ODE45 Sim

[t,w] = ode45(@(t,w) 1/I_r*((tau_s-w*tau_s/w_nl)- Tau_f - param_var.cd*Tau_d*w^2), tspan, w0);

tr_sim = 0;
wterm_sim = 0;

% load varaibles into output variables
w_sim_ar = w;
t_sim_ar = t;

% convert to RPM - not the most efficient code, but easy to read
omega_deg_psec = w*180/pi;
omega_rpm = omega_deg_psec*60/360;
omega_terminal_rpm = omega_rpm(end);  % since we assume the duration is enough to reach terminal velocity, the last point in array is used to define the terminal velocity
    
% Finding terminal velocity from filtered data
wterm_sim = mean(omega_rpm(length(omega_rpm)-20:end)); % [RPM]

% Finding rise time 
threshold = 0.63 * wterm_sim; % rise time = 63% of terminal velocity 
tr_index = find(omega_rpm <= threshold,1,'last'); % index of time value where velocity = 63%
tr_sim = t(tr_index); % rise time value at that index

if pl_flag == true  

        figure(1);
        plot(t,omega_rpm, 'LineWidth',2, 'color', [0 0.4470 0.7410]);
        xline(tr_sim,'r');
        yline(wterm_sim,'r');
        hold on
        
end % end of plot code



