%% Lab Two 
%% Plot One: Positon Data
clear
clc

M1 = readmatrix('EncoderData_NoDelay.txt');

T1 = M1(32:end,1);
C1 = M1(32:end,2);

t1 = T1/(10^6);

figure(1)
plot(t1,C1)
xlabel('Time(sec)')
ylabel('Encoder Counts')
title('Position of Motor Shaft with No Delay')


%% Plot Two: Unfiltered Velocity Data

Y1 = C1/48*2*pi;
for i = 2:length(t1)
    V1(i) = ((Y1(i)- Y1(i-1))/(t1(i)-t1(i-1)));
end

figure(2)
plot(t1,V1*60/(2*pi));
xlabel('Time(sec)');
ylabel('unfiltered velocity (RPM)');
title('Unfiltered Velocity without Delay');


%% Plot Three: Unfiltered Velocity Data with Delay


M2 = readmatrix('EncoderData_10ms.txt');

T2 = M2(18:end,1);
C2 = M2(18:end,2);

Y2 = C2/48*2*pi;

t2 = T2/(10^6);

for i = 2:length(t2)
    V2(i) = ((Y2(i)- Y2(i-1))/(t2(i)-t2(i-1)));
end

figure(3)
plot(t2,V2*60/(2*pi));
xlabel('Time(sec)')
ylabel('unfiltered velocity (RPM)')
title('Unfiltered Velocity with Delay')

%% Plot Four: Filtered Velocity
Tf = [.002,.005,.05]; %sec

tc1 = diff(t1);

% for i = 2:length(V1)
%     delt = diff(t1);
%     %delx = Y1(i)-Y1(i-1);
%     Vf11(i) = (((Y1(i)-Y1(i-1)) + Tf(1)*Vf11(i-1))/(Tf(1)+delt(i-1)))*(60/(2*pi));
% end

Vf11 = 0;
Vf12 = 0;
Vf13 = 0;
for i = 2:length(V1)
    Vf11(i) = (((Y1(i)-Y1(i-1)) + (Tf(1)*Vf11(i-1)))/((Tf(1) + tc1(i-1))));
    Vf12(i) = (((Y1(i)-Y1(i-1)) + (Tf(2)*Vf12(i-1)))/((Tf(2) + tc1(i-1))));
    Vf13(i) = (((Y1(i)-Y1(i-1)) + (Tf(3)*Vf13(i-1)))/((Tf(3) + tc1(i-1))));
end


figure(4)
plot(t1,Vf11*60/(2*pi), Color = [.9 0 .1]);
hold on
plot(t1,Vf12*60/(2*pi), Color = [0 0 1]);
plot(t1,Vf13*60/(2*pi), Color = [0 .8 .2]);
axis([0 2 0 9e3]);
xlabel('Time(sec)')
ylabel('Filtered velocity (RPM)')
title('Filtered Velocity with Delay')
legend('time constant (tc) = .002','tc =.005','tc = .05',Location='southeast')
hold off

%% Plot Five: Config A and B

M3 = readmatrix('EncoderData_NoGearBox.txt');

t3 = M3(21:end, 1)/1e6;
Y3 = M3(21:end, 2)/48*2*pi;

tc3 = diff(t3);

Vf3 = 0;
for i = 2:length(Y3)
    Vf3(i) = (((Y3(i)-Y3(i-1)) + (Tf(3)*Vf3(i-1)))/((Tf(3) + tc3(i-1))));
end

M4 = readmatrix('EncoderData_ConfigB.txt');

t4 = M4(42:end, 1)/1e6;
Y4 = M4(42:end, 2)/48*2*pi;

tc4 = diff(t4);

Vf4 = 0;
for i = 2:length(Y4)
    Vf4(i) = (((Y4(i)-Y4(i-1)) + (Tf(3)*Vf4(i-1)))/((Tf(3) + tc4(i-1))));
end

figure(5)
plot(t3, Vf3*60/(2*pi), Color = [0 0 1]);
hold on
plot(t4, Vf4*60/(2*pi), Color = [.9 0 .1]);
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
title('Config A vs. Config B')
legend('Config A', 'Config B', Location= 'southeast');


%% Plot Five: Config B and C

M5 = readmatrix('EncoderData_ConfigC.txt');

t5 = M5(20:end, 1)/1e6;
Y5 = M5(20:end, 2)/48*2*pi;

tc5 = diff(t5);

Vf5 = 0;
for i = 2:length(Y5)
    Vf5(i) = (((Y5(i)-Y5(i-1)) + (Tf(3)*Vf5(i-1)))/((Tf(3) + tc5(i-1))));
end

figure(6)
plot(t4, Vf4*60/(2*pi), Color = [.9 0 .1]);
hold on
plot(t5, Vf5*60/(2*pi), Color = [0 0 1]);
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
title('Config B vs. Config C')
legend('Config B', 'Config C', Location= 'southeast');

%% Model
%Parameters

% Motor
Jm = 5e-7; % Kgm^2
W_nl = 8200*2*pi/60; % rad/sec
tau_s = 0.0167; % Nm

% Flywheel
d = 4.5*.0254; % m
th = .18*.0254; % m
rho = 1.18*1000; % kg/m^3
v = pi*(d/2)^2*th; % m^3
r = 2*.0254; % m

% nuts and bolts
m = 14/1000; % kg
r2 = .25/2*.0254; % m

% Closed Form Equation for Velocity Profile
tc1 = W_nl*Jm/tau_s;
W = @(t) W_nl*(1-exp(-t/tc1))*60/(2*pi);

% configuration B (Flywheel inertia directly added to motor rotor inertia)
mc = rho*v;
If = 1/2*mc*(d/2)^2;
Ib_a = 1/2*m*r2^2;
Ib_z = m*r^2;
Ib = 8*(Ib_a + Ib_z);
JB = Jm + If + Ib;
tcb = W_nl*JB/tau_s;
WB = @(t) (W_nl*(1-exp((-t/tcb))))*((60)/(2*pi));


%% Plot Six: Config A Experimental and Theoretical

figure(7)
fplot(W,[0 2], Color = [1 0 0]);
hold on
plot(t3, Vf3*60/(2*pi), Color = [0 0 1]);
axis([0 2 0 9000]);
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
title('Config A: Experimental vs. Theoretical');
legend('Theoretical', 'Experimental', Location = 'southeast');

%% Plot Six: Config B Experimental and Theoretical

figure(8)
fplot(WB, [0 10], Color = [1 0 0]);
hold on
plot(t4, Vf4*60/(2*pi), Color = [0 0 1]);
axis([0 10 0 9000]);
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
title('Config B: Experimental vs. Theoretical');
legend('Theoretical', 'Experimental', Location = 'southeast');

%% Personal Plot
M6 = readmatrix('EncoderData_ConfigBSpaced.txt');

t6 = M6(6:end, 1)/1e6;
Y6 = M6(6:end, 2)/48*2*pi;

tc6 = diff(t6);

Vf6 = 0;
for i = 2:length(Y6)
    Vf6(i) = (((Y6(i)-Y6(i-1)) + (Tf(3)*Vf6(i-1)))/((Tf(3) + tc6(i-1))));
end

figure(9)
plot(t4, Vf4*60/(2*pi), Color = [0 0 1]);
hold on
plot(t6, Vf6*60/(2*pi), Color = [.9 0 .1]);
title('Config B vs. Config B Spaced');
xlabel('Time (sec)');
ylabel('Filtered Velocity (RPM)');
legend('Config B', 'Config B Spaced', Location = 'southeast');