%% Motor Simulation with Numerical ODE Approach
%% Analysis of Motor Equation
% The motor equation, also known as the angular acceleration, can be
% derived the torque tau. The two equations of tau used to find the ODE are tau = I_r * a , inertia * angular
% acceleration, and tau = tau_s - w * tau_s/w_nl, the stall torque
% subtracted from  angular velocity * (stall torque /no load speed). By
% setting tau equal to each other, the equation becomes an ode in terms of
% angular acceleration and velocity. By finding the solution to this
% equation, we are given the equation for the angular velocity.

% The inertia term is present within the angular velocity equation. As
% inertia increases, the speed decreases.
%% Simulation

% Simulation parameters

w_nl = 8200*((60)/(2*pi)/10); % (rad/msec) no load speed
tau_s = .17*.0980665; % (Nm) stall torque
I_r = 5*10^-7;  % (Kgm^2) Rotor Intertia
% a = Angular Acceleration
% w = Angular Velocity
% K: Motor Torque Constant
% R: Electrical resistance 

% Motor Equation

    % tau = (1/R) * ((K * V) - (K^2 * w_nl))
    % tau = tau_s - w * tau_s/w_nl
    % tau = I_r * a
       
    % combination of the two equations above give and ODE
    % a = (1/I_r) * (tau_s - w * (tau_s/w_nl))

%% ODE Method

tspan = [0 15]; % span of time
w0 = 0; % initial condition for angular velocity

% using ode45 command to find solution to angular acceleration ode

[t,w] = ode45(@(t,w) ((1/I_r) * (tau_s - w * (tau_s/w_nl)))*((2*pi/(60))), tspan, w0);

1/I_r * ((taus_s - w * tau_s/w_nl) - 

figure(1)
plot(t,w)
title('Motor Speed')
xlabel('time (sec)')
ylabel('Angular Velocity (RPM)')

