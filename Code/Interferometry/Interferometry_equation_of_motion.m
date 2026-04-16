clear;
close all;

%% Input parameters

d = 5000e3;          % Initial distance (m)
dt = 0.01;           % Time step
time_span = 1e5;     % Time span of the simulation (s)
G = 6.67384e-11;     % Gravitational constant
M = 35342917353;     % Asteroid mass [kg]
v0 = 100;            % x'(0) [m/s]
h = 1.15*1000;       % Altitude [m]
S = 220000;          % Satellite separation [m]


%% Function definition
mu = G*M;

ode_system = @(t, y) [y(3);                       % dx/dt
                      y(4);                       % dy/dt
                      -mu/(norm(y(1:3))^3)*y(1);  % dvx/dt
                      -mu/(norm(y(1:3))^3)*y(2)]; % dvy/dt


%% Initial conditions
y01 =[-d h v0 0]';
y02 =[-d+S h v0 0]';


%% Solver conf.
% tspan = [0 time_span];
tspan = 0:dt:time_span;
options = odeset( ...
    'RelTol',1e-14, ...
    'AbsTol',1e-14);


%% Solver
[t1,y1] = ode45(ode_system, tspan, y01, options);
[t2,y2] = ode45(ode_system, tspan, y02, options);


%% Extract solutions
x1 = y1(:,1);
v1 = y1(:,3);
x2 = y2(:,1);
v2 = y2(:,3);

diffv = v2-v1;

[max_val, i_max] = max(diffv);
[min_val, i_min] = min(diffv);

semi_T = t1(i_min) - t1(i_max);

f = 1/(2*semi_T);

disp('Frequency [mHz]:')
disp(f*1000);

disp('Max diff:')
disp(max_val);

%% Plot
figure
plot(y1(:,1),y1(:,2));
hold on
scatter(0,0,'o')

figure
plot(t1,x1,'LineWidth',2)
hold on
plot(t1,x2,'LineWidth',2)
xlabel('t [s]')
ylabel('x(t) [m]')
title('Satellite positions')
legend('Sat 1', 'Sat 2')
grid on

figure
plot(t1,v1,'LineWidth',1)
hold on
plot(t2,v2,'LineWidth',1)
xlabel('t [s]')
ylabel('v [m/s]')
title('Velocity')
legend('Sat 1', 'Sat 2')
grid on

figure
plot(t1,diffv,'LineWidth',1)
xlabel('t [s]')
ylabel('Separation rate [m/s]')
title('Separation rate')
legend('Sat 1', 'Sat 2')
grid on