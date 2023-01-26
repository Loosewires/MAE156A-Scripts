%% Lab 3: Plots
close all
clear
clc

%% Parameters
    J_m = 5e-7; % Rotational inertia of motor shaft [kgm^2]
    W_nl = 8200*2*pi/60; % No load speed [rad/sec]
    Tau_s = 0.0167; % Stall torque [Nm]
    K = .014; % Motor constant
    R = .014/.0014; % Motor resistance
    Tf = .05; % Filtering constant

%% Config 1

    % 1.1 - 50% PWM
        % Closed Form Solution
    
        Tau_s_11 = K/R*(12*.5); % Stall torque [Nm]
        W_nl_11 = 12*.5/K; % No load speed [rad/s]
        J_11 = 7.34e-6; % Config 1 rotational inertia [kgm^2]
        Tau_fm_11 = 4.69e-4; % Motor torque required to overcome friction [Nm]
        W_nl_11 = W_nl_11*(Tau_s_11-Tau_fm_11)/Tau_s_11; % No load speed with friction [rad/s]
        TC_11 = W_nl_11 * J_11 / Tau_s_11; % Config 1 time constant
        W_11 = @(t) W_nl_11*(1-exp(-t/TC_11))*60/(2*pi); % Closed form solution

        % Experimental Data

        M11 = readmatrix('EncoderData_Config11.txt');

        t11 = M11(1:end, 1)/1e6;
        Y11 = M11(1:end, 2)/48*2*pi;
        
        tc11 = diff(t11);
        
        Vf11 = 0;
        for i = 2:length(Y11)
            Vf11(i) = (((Y11(i)-Y11(i-1)) + (Tf*Vf11(i-1)))/((Tf + tc11(i-1))));
        end

        % Plot
        figure(1)
        fplot(W_11, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t11, Vf11*60/(2*pi), Color = [0 0 1]);
        title('Config 1.1: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity [Tf = 0.05] (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 4500]);

    % 1.2 - 75% PWM
        % Closed Form Solution

        Tau_s_12 = K/R*(12*.75); % Stall torque [Nm]
        W_nl_12 = 12*.75/K; % No load speed [rad/s]
        J_12 = J_11; % Config 1 rotational inertia [kgm^2]
        Tau_fm_12 = Tau_fm_11; % Motor torque required to overcome friction [Nm]
        W_nl_12 = W_nl_12*(Tau_s_12-Tau_fm_12)/Tau_s_12; % No load speed with friction [rad/s]
        TC_12 = W_nl_12 * J_12 / Tau_s_12; % Config 1 time constant
        W_12 = @(t) W_nl_12*(1-exp(-t/TC_12))*60/(2*pi); % Closed form solution

        % Experimental Data

        M12 = readmatrix('EncoderData_Config12.txt');

        t12 = M12(1:end, 1)/1e6;
        Y12 = M12(1:end, 2)/48*2*pi;
        
        tc12 = diff(t12);
        
        Vf12 = 0;
        for i = 2:length(Y12)
            Vf12(i) = (((Y12(i)-Y12(i-1)) + (Tf*Vf12(i-1)))/((Tf + tc12(i-1))));
        end


        % Plot
        figure(2)
        fplot(W_12, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t12, Vf12*60/(2*pi), Color = [0 0 1]);
        title('Config 1.2: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 6300]);

    % 1.3 - 100% PWM
        % Closed Form Solution

        Tau_s_13 = K/R*(12*1); % Stall torque [Nm]
        W_nl_13 = 12*1/K; % No load speed [rad/s]
        J_13 = J_11; % Config 1 rotational inertia [kgm^2]
        Tau_fm_13 = Tau_fm_11; % Motor torque required to overcome friction [Nm]
        W_nl_13 = W_nl_13*(Tau_s_13-Tau_fm_13)/Tau_s_13; % No load speed with friction [rad/s]
        TC_13 = W_nl_13 * J_13 / Tau_s_13; % Config 1 time constant
        W_13 = @(t) W_nl_13*(1-exp(-t/TC_13))*60/(2*pi); % Closed form solution

        % Experimental Data

        M13 = readmatrix('EncoderData_Config13.txt');

        t13 = M13(1:end, 1)/1e6;
        Y13 = M13(1:end, 2)/48*2*pi;
        
        tc13 = diff(t13);
        
        Vf13 = 0;
        for i = 2:length(Y13)
            Vf13(i) = (((Y13(i)-Y13(i-1)) + (Tf*Vf13(i-1)))/((Tf + tc13(i-1))));
        end

        % Plot
        figure(3)
        fplot(W_13, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t13, Vf13*60/(2*pi), Color = [0 0 1]);
        title('Config 1.3: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 8200]);

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

        % Experimental Data

        M21 = readmatrix('EncoderData_Config21.txt');

        t21 = M21(1:end, 1)/1e6;
        Y21 = M21(1:end, 2)/48*2*pi;
        
        tc21 = diff(t21);
        
        Vf21 = 0;
        for i = 2:length(Y21)
            Vf21(i) = (((Y21(i)-Y21(i-1)) + (Tf*Vf21(i-1)))/((Tf + tc21(i-1))));
        end


        % Plot
        figure(4)
        fplot(W_21, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t21, Vf21*60/(2*pi), Color = [0 0 1]);
        title('Config 2.1: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 4500]);

    % 2.2 - 75% PWM
        % Closed Form Solution

        Tau_s_22 = K/R*(12*.75); % Stall torque [Nm]
        W_nl_22 = 12*.75/K; % No load speed [rad/s]
        J_22 = J_21; % Config 1 rotational inertia [kgm^2]
        Tau_fm_22 = Tau_fm_21; % Motor torque required to overcome friction [Nm]
        W_nl_22 = W_nl_22*(Tau_s_22-Tau_fm_22)/Tau_s_22; % No load speed with friction [rad/s]
        TC_22 = W_nl_22 * J_22 / Tau_s_22; % Config 1 time constant
        W_22 = @(t) W_nl_22*(1-exp(-t/TC_22))*60/(2*pi); % Closed form solution

        % Experimental Data

        M22 = readmatrix('EncoderData_Config22.txt');

        t22 = M22(1:end, 1)/1e6;
        Y22 = M22(1:end, 2)/48*2*pi;
        
        tc22 = diff(t22);
        
        Vf22 = 0;
        for i = 2:length(Y22)
            Vf22(i) = (((Y22(i)-Y22(i-1)) + (Tf*Vf22(i-1)))/((Tf + tc22(i-1))));
        end

        % Plot
        figure(5)
        fplot(W_22, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t22, Vf22*60/(2*pi), Color = [0 0 1]);
        title('Config 2.2: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity [Tf = 0.05] (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 6200]);

    % 2.3 - 100% PWM
        % Closed Form Solution

        Tau_s_23 = K/R*(12*1); % Stall torque [Nm]
        W_nl_23 = 12*1/K; % No load speed [rad/s]
        J_23 = J_21; % Config 1 rotational inertia [kgm^2]
        Tau_fm_23 = Tau_fm_21; % Motor torque required to overcome friction [Nm]
        W_nl_23 = W_nl_23*(Tau_s_23-Tau_fm_23)/Tau_s_23; % No load speed with friction [rad/s]
        TC_23 = W_nl_23 * J_23 / Tau_s_23; % Config 1 time constant
        W_23 = @(t) W_nl_23*(1-exp(-t/TC_23))*60/(2*pi); % Closed form solution

        % Experimental Data

        M23 = readmatrix('EncoderData_Config23.txt');

        t23 = M23(1:end, 1)/1e6;
        Y23 = M23(1:end, 2)/48*2*pi;
        
        tc23 = diff(t23);
        
        Vf23 = 0;
        for i = 2:length(Y23)
            Vf23(i) = (((Y23(i)-Y23(i-1)) + (Tf*Vf23(i-1)))/((Tf + tc23(i-1))));
        end

        % Plot
        figure(6)
        fplot(W_23, [0 5], Color = [.9 0 .1]);
        hold on
        plot(t23, Vf23*60/(2*pi), Color = [0 0 1]);
        title('Config 2.3: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 5 0 8300]);

%% Config 3

    % 3.1 - 50% PWM
        % Closed Form Solution

        Tau_s_31 = K/R*(12*.5); % Stall torque [Nm]
        W_nl_31 = 12*.5/K; % No load speed [rad/s]
        J_31 = 1.95e-5; % Config 3.1 rotational inertia [kgm^2]
        Tau_fm_31 = 9.89e-4; % Motor torque required to overcome friction [Nm]
        W_nl_31 = W_nl_31*(Tau_s_31-Tau_fm_31)/Tau_s_31; % No load speed with friction [rad/s]
        TC_31 = W_nl_31 * J_31 / Tau_s_31; % Config 3.1 time constant
        W_31 = @(t) W_nl_31*(1-exp(-t/TC_31))*60/(2*pi); % Closed form solution

        % Experimental Data

        M31 = readmatrix('EncoderData_Config31.txt');

        t31 = M31(1:end, 1)/1e6;
        Y31 = M31(1:end, 2)/48*2*pi;
        
        tc31 = diff(t31);
        
        Vf31 = 0;
        for i = 2:length(Y31)
            Vf31(i) = (((Y31(i)-Y31(i-1)) + (Tf*Vf31(i-1)))/((Tf + tc31(i-1))));
        end

        % Plot
        figure(7)
        fplot(W_31, [0 10], Color = [.9 0 .1]);
        hold on
        plot(t31, Vf31*60/(2*pi), Color = [0 0 1]);
        title('Config 3.1: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity [Tf = 0.05] (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 10 0 4500]);

    % 3.2 - 75% PWM
        % Closed Form Solution

        Tau_s_32 = K/R*(12*.75); % Stall torque [Nm]
        W_nl_32 = 12*.75/K; % No load speed [rad/s]
        J_32 = J_31; % Config 3.3 rotational inertia [kgm^2]
        Tau_fm_32 = Tau_fm_31; % Motor torque required to overcome friction [Nm]
        W_nl_32 = W_nl_32*(Tau_s_32-Tau_fm_32)/Tau_s_32; % No load speed with friction [rad/s]
        TC_32 = W_nl_32 * J_32 / Tau_s_32; % Config 3.3 time constant
        W_32 = @(t) W_nl_32*(1-exp(-t/TC_32))*60/(2*pi); % Closed form solution

        % Experimental Data

        M32 = readmatrix('EncoderData_Config32.txt');

        t32 = M32(1:end, 1)/1e6;
        Y32 = M32(1:end, 2)/48*2*pi;
        
        tc32 = diff(t32);
        
        Vf32 = 0;
        for i = 2:length(Y32)
            Vf32(i) = (((Y32(i)-Y32(i-1)) + (Tf*Vf32(i-1)))/((Tf + tc32(i-1))));
        end

        % Plot
        figure(8)
        fplot(W_32, [0 15], Color = [.9 0 .1]);
        hold on
        plot(t32, Vf32*60/(2*pi), Color = [0 0 1]);
        title('Config 3.2: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 10 0 6200]);

    % 3.3 - 100% PWM
        % Closed Form Solution

        Tau_s_33 = K/R*(12*1); % Stall torque [Nm]
        W_nl_33 = 12*1/K; % No load speed [rad/s]
        J_33 = J_31; % Config 3.3 rotational inertia [kgm^2]
        Tau_fm_33 = Tau_fm_31; % Motor torque required to overcome friction [Nm]
        W_nl_33 = W_nl_33*(Tau_s_33-Tau_fm_33)/Tau_s_33; % No load speed with friction [rad/s]
        TC_33 = W_nl_33 * J_33 / Tau_s_33; % Config 3.3 time constant
        W_33 = @(t) W_nl_33*(1-exp(-t/TC_33))*60/(2*pi); % Closed form solution

        % Experimental Data

        M33 = readmatrix('EncoderData_Config33.txt');

        t33 = M33(1:end, 1)/1e6;
        Y33 = M33(1:end, 2)/48*2*pi;
        
        tc33 = diff(t33);
        
        Vf33 = 0;
        for i = 2:length(Y33)
            Vf33(i) = (((Y33(i)-Y33(i-1)) + (Tf*Vf33(i-1)))/((Tf + tc33(i-1))));
        end

        % Plot
        figure(9)
        fplot(W_33, [0 15], Color = [.9 0 .1]);
        hold on
        plot(t33, Vf33*60/(2*pi), Color = [0 0 1]);
        title('Config 3.3: Theoretical vs. Experimental');
        xlabel('Time (sec)');
        ylabel('Filtered Velocity [Tf = 0.05] (RPM)');
        legend('Theoretical', 'Experimental', Location = 'southeast');
        axis([0 10 0 8200]);



            