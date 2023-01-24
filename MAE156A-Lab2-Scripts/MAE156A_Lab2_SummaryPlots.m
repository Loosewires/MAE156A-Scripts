%% Lab 2: Summary Plots

%% Parameters

    % Motor
        J_m = 5e-7; % Kgm^2
        W_nl = 8200*2*pi/60; % rad/sec
        Tau_s = 0.0167; % Nm

    % Config A
        J_A = 5e-7; % Nm
        TC_A = W_nl * J_A / Tau_s; % s
        W_A = @(t) W_nl*(1-exp(-t/TC_A))*60/(2*pi);

    % Config B
        J_B = 1.95e-5; % Nm
        TC_B = W_nl * J_B / Tau_s; % s
        W_B = @(t) W_nl*(1-exp(-t/TC_B))*60/(2*pi);

    % Config C
        J_C = 1.88e-5; % Nm
        TC_C = W_nl * J_C / Tau_s; % s
        W_C = @(t) W_nl*(1-exp(-t/TC_C))*60/(2*pi);
    
    % Time Constant
        Tf = .05; % s

%% Plot 1: Config A Experimental and Theoretical

M3 = readmatrix('EncoderData_NoGearBox.txt');

t3 = M3(21:end, 1)/1e6;
Y3 = M3(21:end, 2)/48*2*pi;

tc3 = diff(t3);

Vf3 = 0;
for i = 2:length(Y3)
    Vf3(i) = (((Y3(i)-Y3(i-1)) + (Tf*Vf3(i-1)))/((Tf + tc3(i-1))));
end

figure(1)
fplot(W_A, [0 2], Color = [.9 0 .1]);
hold on
plot(t3, Vf3*60/(2*pi), Color = [0 0 1]);
title('Config A: Experimental vs Theoretical');
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
legend('Theoretical', 'Experimental', Location = 'southeast');
axis([0 2 0 9000]);

CA_TermVel_Th = sum(W_A(1.9:.001:2))/length(W_A(1.9:.001:2));
CA_TermVel_Ex = sum(Vf3(1679:1776))*60/(2*pi)/(1776-1679);

CA_RiseTime_Th = -TC_A*log(1-.63*CA_TermVel_Th*2*pi/(W_nl*60));
CA_RiseTime_Ex = t3(113);

%% Plot 2: Config B Experimental and Theoretical

M4 = readmatrix('EncoderData_ConfigB.txt');

t4 = M4(42:end, 1)/1e6;
Y4 = M4(42:end, 2)/48*2*pi;

tc4 = diff(t4);

Vf4 = 0;
for i = 2:length(Y4)
    Vf4(i) = (((Y4(i)-Y4(i-1)) + (Tf*Vf4(i-1)))/((Tf + tc4(i-1))));
end

figure(2)
fplot(W_B, [0 5], Color = [.9 0 .1]);
hold on
plot(t4, Vf4*60/(2*pi), Color = [0 0 1]);
title('Config B: Experimental vs Theoretical');
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
legend('Theoretical', 'Experimental', Location = 'southeast');
axis([0 5 0 9000]);

CB_TermVel_Th = sum(W_B(4.9:.001:5))/length(W_B(4.9:.001:5));
CB_TermVel_Ex = sum(Vf4(4141:4226))*60/(2*pi)/(4226-4141);

CB_RiseTime_Th = -TC_B*log(1-.63*CB_TermVel_Th*2*pi/(W_nl*60));
CB_RiseTime_Ex = t4(955);

%% Plot 3: Config C Experimental and Theoretical

M5 = readmatrix('EncoderData_ConfigC.txt');

t5 = M5(20:end, 1)/1e6;
Y5 = M5(20:end, 2)/48*2*pi;

tc5 = diff(t5);

Vf5 = 0;
for i = 2:length(Y5)
    Vf5(i) = (((Y5(i)-Y5(i-1)) + (Tf*Vf5(i-1)))/((Tf + tc5(i-1))));
end

figure(3)
fplot(W_C, [0 5], Color = [.9 0 .1]);
hold on
plot(t5, Vf5*60/(2*pi), Color = [0 0 1]);
title('Config C: Experimental vs Theoretical');
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
legend('Theoretical', 'Experimental', Location = 'southeast');
axis([0 5 0 9000]);

CC_TermVel_Th = sum(W_C(4.9:.001:5))/length(W_A(4.9:.001:5));
CC_TermVel_Ex = sum(Vf5(4119:4197))*60/(2*pi)/(4197-4119);

CC_RiseTime_Th = -TC_C*log(1-.63*CC_TermVel_Th*2*pi/(W_nl*60));
CC_RiseTime_Ex = t5(860);