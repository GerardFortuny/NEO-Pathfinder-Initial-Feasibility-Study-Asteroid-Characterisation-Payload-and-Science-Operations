clear all
close all

%% Data
% SU49
m = 7.3E+10;          % Asteroid mass [kg]
d = 600;              % Asteroid diameter [m]

% % SV19
% m = 1E+8;               % Asteroid mass [kg]
% d = 43;

% SC data
sigma_v0 = 3.14E-04;  % Doppler velocity error [m/s] for 1 s integration in X-Band
t_int = 600;          % Integration time [s]
G = 6.6743E-11;       % Gravitational constant [m^3 kg^-1 s^-2]



mu = G*m;

h_vec = 0:10:10000;    % Altitude from surface [m]
v_vec = 0:0.1:5;       % Flyby velocity (periapsis) [m/s]


%% Calc

sigma_v = sigma_v0 / sqrt(t_int);

for i = 1:length(h_vec)
    centr_dist = d/2 + h_vec(i);

    for j = 1:length(v_vec)

        mass_err(i,j) = sigma_v * (centr_dist * v_vec(j) / (mu)) * 100;

        % Velocity resolution check
        v_inf = sqrt(v_vec(j)^2 - 2*mu/centr_dist);
        delta_V_t(i,j) = abs(v_inf-v_vec(j));

    end

    v_escape(i) = sqrt(2*G*m/centr_dist);
end


%% Plot

[X, Y] = meshgrid(v_vec, h_vec);

figure
p = pcolor(X, Y, mass_err);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
shading interp
hold on
contour(X, Y, mass_err, [2.5 2.5], 'r', 'LineWidth', 2, 'DisplayName','Mass error = 2.5%')
hold on
plot(v_escape,h_vec,'g','LineWidth', 2, 'DisplayName','Escape velocity')
colorbar
xlabel('Flyby velocity [m/s]')
ylabel('Altitude [m]')
title('Mass error [%]')
legend show


figure
p = pcolor(X, Y, delta_V_t);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
shading interp
hold on
contour(X, Y, delta_V_t, [(0.32e-3)/sqrt(t_int) (0.32e-3)/sqrt(t_int)], 'r', 'LineWidth', 2, 'DisplayName','Velocity difference = velocity error')
colorbar
xlabel('Flyby velocity [m/s]')
ylabel('Altitude [m]')
title('Mass error [%]')
legend show

