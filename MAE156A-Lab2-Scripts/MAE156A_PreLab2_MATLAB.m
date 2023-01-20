%% parameters

    % motor
    Jm = 5e-7; % Kgm^2
    W_nl = 8200*2*pi/60; % rad/sec
    tau_s = 0.0167; % Nm
        
    % Flywheel
    d = 4.5*.0254; % m
    th = .18*.0254/100; % m
    rho = 1.18*1000; % kg/m^3
    v = pi*(d/2)^2*th; %m^3
    r = 2*.0254; % m
    
    % nuts and bolts
    m = 14/1000; % kg
    r2 = .25/2*.0254; % m

%% Closed Form Equation for V elocity Profile
tc = W_nl*Jm/tau_s*1000;
W = @(t) W_nl*(1-exp(-t/tc))*60/(2*pi);

    % configuration A (No inertia added to motor rotor)
    figure(1)
    fplot(W,[0, 500]);
    xlabel("Time (milliseconds)")
    ylabel("Motor Velocity (RPM)")
    title("Configuration A")

    % configuration B (Flywheel inertia directly added to motor rotor inertia)
    mc = rho*v;
    If = 1/2*mc*(d/2)^2;
    Ib_a = 1/2*m*r2^2;
    Ib_z = m*r^2;
    Ib = 8*(Ib_a + Ib_z);
    JB = Jm + If + Ib;
    tcb = W_nl*JB/tau_s*1000;
    WB = @(t) (W_nl*(1-exp((-t/tcb))))*((60)/(2*pi));
    
    figure(2)
    fplot(WB,[0, 100000]);
    xlabel("Time (milliseconds)");
    ylabel("Motor Velocity (RPM)");
    title("Configuration B");
