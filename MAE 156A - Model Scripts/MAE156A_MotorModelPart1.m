%% Lab 3: Plots

%% Parameters
    J_m = 5e-7; % Rotational inertia of motor shaft [kgm^2]
    W_nl = 8200*2*pi/60; % No load speed [rad/sec]
    Tau_s = 0.0167; % Stall torque [Nm]
    K = .014; % Motor constant
    R = .014/.0014; % Motor resistance

%% Config 1

    % 1.1 - 50% PWM
        % Closed Form Solution
    
        Tau_s_11 = K/R*6; % Stall torque [Nm]
        W_nl_11 = 6/K; % No load speed [rad/s]
        J_11 = 7.34e-6; % Config 1 rotational inertia [kgm^2]
        Tau_fm_11 = 4.69e-4; % Motor torque required to overcome friction [Nm]
        W_nl_11 = W_nl_11*(Tau_s_11-Tau_fm_11)/Tau_s_11; % No load speed with friction [rad/s]
        TC_11 = W_nl_11 * J_11 / Tau_s_11; % Config 1 time constant
        W_11 = @(t) W_nl_11*(1-exp(-t/TC_11))*60/(2*pi); % Closed form solution

        % Plot
        figure(1)
        fplot(W_11, [0 5]);

    % 1.2 - 75% PWM
        % Closed Form Solution

        Tau_s_12 = K/R*(12*.75); % Stall torque [Nm]
        W_nl_12 = 12*.75/K; % No load speed [rad/s]
        J_12 = J_11; % Config 1 rotational inertia [kgm^2]
        Tau_fm_12 = Tau_fm_11; % Motor torque required to overcome friction [Nm]
        W_nl_12 = W_nl_12*(Tau_s_12-Tau_fm_12)/Tau_s_12; % No load speed with friction [rad/s]
        TC_12 = W_nl_12 * J_12 / Tau_s_12; % Config 1 time constant
        W_12 = @(t) W_nl_12*(1-exp(-t/TC_12))*60/(2*pi); % Closed form solution

        % Plot
        figure(2)
        fplot(W_12, [0 5]);

    % 1.3 - 100% PWM
        % Closed Form Solution

        Tau_s_13 = K/R*(12*1); % Stall torque [Nm]
        W_nl_13 = 12*1/K; % No load speed [rad/s]
        J_13 = J_11; % Config 1 rotational inertia [kgm^2]
        Tau_fm_13 = Tau_fm_11; % Motor torque required to overcome friction [Nm]
        W_nl_13 = W_nl_13*(Tau_s_13-Tau_fm_13)/Tau_s_13; % No load speed with friction [rad/s]
        TC_13 = W_nl_13 * J_13 / Tau_s_13; % Config 1 time constant
        W_13 = @(t) W_nl_13*(1-exp(-t/TC_13))*60/(2*pi); % Closed form solution

        % Plot
        figure(3)
        fplot(W_13, [0 5]);

%% Config 2

    % 2.1 - 50% PWM
        % Closed Form Solution
        Tau_s_21 = K/R*(12*.5); % Stall torque [Nm]
        W_nl_21 = 12*.5/K; % No load speed [rad/s]
        J_21 = 1.03e-5; % Config 1 rotational inertia [kgm^2]
        Tau_fm_21 = 5.94e-4; % Motor torque required to overcome friction [Nm]
        W_nl_21 = W_nl_21*(Tau_s_21-Tau_fm_21)/Tau_s_21; % No load speed with friction [rad/s]
        TC_21 = W_nl_21 * J_21 / Tau_s_21; % Config 1 time constant
        W_21 = @(t) W_nl_21*(1-exp(-t/TC_21))*60/(2*pi); % Closed form solution

        % Plot
        figure(4)
        fplot(W_21, [0 5]);

    % 2.2 - 75% PWM
        % Closed Form Solution

        Tau_s_22 = K/R*(12*.75); % Stall torque [Nm]
        W_nl_22 = 12*.75/K; % No load speed [rad/s]
        J_22 = J_21; % Config 1 rotational inertia [kgm^2]
        Tau_fm_22 = Tau_fm_21; % Motor torque required to overcome friction [Nm]
        W_nl_22 = W_nl_22*(Tau_s_22-Tau_fm_22)/Tau_s_22; % No load speed with friction [rad/s]
        TC_22 = W_nl_22 * J_22 / Tau_s_22; % Config 1 time constant
        W_22 = @(t) W_nl_22*(1-exp(-t/TC_22))*60/(2*pi); % Closed form solution

        % Plot
        figure(5)
        fplot(W_22, [0 5]);

    % 2.3 - 100% PWM
        % Closed Form Solution

        Tau_s_23 = K/R*(12*1); % Stall torque [Nm]
        W_nl_23 = 12*1/K; % No load speed [rad/s]
        J_23 = J_11; % Config 1 rotational inertia [kgm^2]
        Tau_fm_23 = Tau_fm_21; % Motor torque required to overcome friction [Nm]
        W_nl_23 = W_nl_23*(Tau_s_23-Tau_fm_23)/Tau_s_23; % No load speed with friction [rad/s]
        TC_23 = W_nl_23 * J_23 / Tau_s_23; % Config 1 time constant
        W_23 = @(t) W_nl_23*(1-exp(-t/TC_23))*60/(2*pi); % Closed form solution

        % Plot
        figure(6)
        fplot(W_23, [0 5]);

%% Config 3

    % 3.1 - 50% PWM


            