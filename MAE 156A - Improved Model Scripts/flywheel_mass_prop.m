function [J_eff, mfw, Tau_f] = flywheel_mass_prop(config, param_var, param_fixed)
%
% This function calcualte the mass properties of the flywheel and uses the value to return
% the effective inertia of the system and the mass of the flywheel.
% Assumption: The nuts and bolts can be treated as point masses.

% Parameters
rs = .002; % Radius of output shaft [m]
g = 9.8; % Acceleration due to gravity [m/s^2]
Lc = .0099; % Length from load to left part [m]
Lb = .0066; % length from left part to right part [m]

nbolt = 0;
for i = 1:length(config.nut_ar)
    if config.nut_ar(i) == 0
        
    else
        nbolt = nbolt + 1;
    end
end

nnut = sum(config.nut_ar); % number of nuts

mfw_nobolts = param_fixed.fw_mass; % [kg]
    
J_disk = 7.58e-5; % inertia of flywheel disk [kg*m^2]
J_bolt = 2e-5; % inertia due to bolts [kg*m^2]
J_hub = 1.73e-7; % inertia due to hub [kg*m^2]
J_nut = 8.26e-6; % inertia due to nuts [kg*m^2]
J_motor_internal = param_var.jm; % inertia of internal motor rotor [kg*m^2]

% Total mass of flywheel + all components
mfw = mfw_nobolts + nbolt*param_fixed.bolt_mass + nnut*param_fixed.nut_mass + param_fixed.hub_mass;

% Total inertia of flywheel + all components in [kg*m^2]
J_eff = (J_disk + J_hub + nbolt*J_bolt + nnut*J_nut) / ((param_fixed.ngear^2)) + J_motor_internal;

% Calculate Friction
Tau_f = param_var.mu*rs*g*mfw/param_fixed.ngear * (2*(Lc+Lb)/Lb - 1);

end