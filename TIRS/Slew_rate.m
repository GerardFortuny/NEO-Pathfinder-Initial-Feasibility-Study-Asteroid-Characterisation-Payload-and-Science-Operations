clear all
close all


%% General data
FOV = 6.35e-3; % [rad]
R_a = 300; % Asteroid radius [m]
t_int = 1.8; % OTES integration time [s]



%% Calculations first part

j = 1;
for rot_period = 2:0.25:10
    i = 1;
    for h = 0:10:5000

        % Lambda of puntual coverture
        lambd = - FOV/2 + asin((R_a+h)*sin(FOV/2)/R_a);

        if (R_a + h) * sin(FOV/2) / R_a > 1
            error('Error.');
        end

        % Time for one pass:
        rot = 2*pi/(rot_period*3600); % Asteroid rotation [rad/s]
        t_p = 2*lambd/rot;

        % Angle of half pass (taking into account the margin for calibration)
        lambd_max = acos(R_a/(R_a+h));
        S_max = pi/2 + FOV/2 - lambd_max;

        % Slew rate
        SR(i,j) = rad2deg(2*S_max/t_p); %[rad/s]

        % Calculate surface coverage
        Cov(i) = sin(acos(R_a/(R_a+h)))*100;

        h_plot(i) = h;

        i = i+1;


    end
    disp(rot_period)

    rot_plot(j) = rot_period;

    j = j+1;
end

%% Plot

figure
plot(h_plot,Cov,LineWidth=2);
ylabel('Surface coverage [%]')
xlabel('Altitude [m]')
grid on

[X, Y] = meshgrid(rot_plot, h_plot);
figure
hold on
% contour3(X, Y, SR, [1 1], 'r', 'LineWidth', 2);
surf(X, Y, SR)
pcolor(X, Y, SR)
shading interp
xlabel('Rotation period [h]')
ylabel('Altitude [m]')
zlabel('SR')
colorbar
% clim([0 6])
caxis([0 6])
legend('Slew rate = 1 º/s')
title('Slew rate [º/s]')






