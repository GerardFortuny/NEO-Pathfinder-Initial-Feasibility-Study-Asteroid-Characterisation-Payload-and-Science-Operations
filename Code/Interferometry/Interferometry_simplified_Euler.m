clear all
close all


%% Precision configuration

digits(20)   % Significant digits


%% Input parameters

x0 = vpa(-10e6);          % Initial distance (m)
d = 0;
dt = vpa(0.1);            % Time step [s]
G = vpa(6.67384e-11);     % Gravitational constant
M = vpa(35342917353);     % Apophis mass [kg]
v0 = vpa(5800);           % Initial velocity [m/s]
h = vpa(150.15*1000);     % Altitude [m]
S = vpa(220000);          % Satellite separation [m]



%% Initialisation

x1 = x0;
x2 = x0 + S;
t = vpa(0);
v1 = v0;
v2 = v0;


i = 1;

%% Main loop

while x1 < -x0

    % Store values for plot
    x1_list(i) = double(x1 / 1000);
    v1_list(i) = double(v1);
    v2_list(i) = double(v2);
    t_list(i) = double(t);
    x2_list(i) = double(x2 / 1000);
    x_mean(i) = double((x1 + x2) / 2000);
    diff(i) = double((x2 - x1) - S);
    dt_list(i) = double(dt);

    % Compute accelerations

    a_x1 = vpa(G*M/(((d-x1)^2)+h^2)*((d-x1)/sqrt(((d-x1)^2)+h^2)));

    a_x2 = vpa(G*M/(((d-x2)^2)+h^2)*((d-x2)/sqrt(((d-x2)^2)+h^2)));

    % Update velocities and positions
    v1 = v1 + a_x1 * dt;
    v2 = v2 + a_x2 * dt;

    x1 = x1 + v1 * dt;
    x2 = x2 + v2 * dt;
    t = t + dt;

    i = i + 1;

    disp(double(x1 / 1000));
end


% Compute velocity differences


diffv = -v1_list + v2_list;

% Calculate frequency

[max_val, i_max] = max(diffv);
[min_val, i_min] = min(diffv);

semi_T = t_list(i_min) - t_list(i_max);

f = vpa(1/(2*semi_T),5);

%% Plots

disp('Frequency [mHz]:')
disp(f*1000);

disp('Max diff:')
disp(max_val);

figure;
plot(x_mean, diff);
xlabel('x mean [km]');
ylabel('Distance difference [m]');
grid on;
xlim([(double(x0 + S/2) / 1000) -(double(x0 - S/2) / 1000)]);

figure;
plot(x_mean, dt_list);
xlabel('x mean [km]');
ylabel('dt [s]');
grid on;
xlim([(double(x0 + S/2) / 1000) -(double(x0 - S/2) / 1000)]);

figure;
plot(x_mean, diffv, 'LineWidth', 2);
xlabel('x mean [km]');
ylabel('Separation rate [m/s]');
grid on;
xlim([(double(x0 + S/2) / 1000) -(double(x0 - S/2) / 1000)]);


figure;
plot(t_list, diffv, 'LineWidth', 2);
xlabel('t [s]');
ylabel('Separation rate [m/s]');
grid on;
xlim([min(t_list) max(t_list)]);

