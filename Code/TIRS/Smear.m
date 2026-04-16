clear all
close all


%% General data
FOV = 6.35e-3; % [rad]
R_a = 300; % Asteroid radius [m]
t_int = 1.8; % OTES integration time [s]


%% Calculations

j = 1;

for rot_period = 2:0.25:10 % Asteroid rotation period [h]
    i = 1;
    for h = 1000:50:15000
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
        SR = 2*S_max/t_p; %[rad/s]

        % Calculate scaning velocity
        m = 1;
        for d = 0:1:cos(pi/2 - lambd_max)*R_a % We look at all latitudes to see where we have the worst ratio
           
            r = sqrt(R_a^2 - d^2);
            lambd_s = atan2(d,r); %pi/2 - atan(r/d);
            s_s = sqrt(R_a^2 + (R_a + h)^2 -2*R_a*(R_a+h)*cos(lambd_s)); %Cos theorem
            eta_s = asin(R_a*sin(lambd_s)/s_s); %Sin theorem
            del_lambd = lambd_s - asin((R_a+h)*sin(eta_s - FOV)/R_a) + (eta_s - FOV);
            D = del_lambd*R_a;

            lambd_c = asin((R_a+h)*sin(eta_s - FOV/2)/R_a) - (eta_s - FOV/2);
            s_c = sqrt(R_a^2 + (R_a + h)^2 -2*R_a*(R_a+h)*cos(lambd_c));

            v_rot = r*2*pi/(rot_period*3600); % Asteroid surface rotation [m/s]

            v_slew = s_c*SR;

            v_T = sqrt(v_rot^2 + v_slew^2);

            ratio(m) = 100*v_T*t_int/D;

            m = m+1;

        end

        ratio_v(i,j) = max(ratio);

        h_plot(i) = h;

        i = i+1;

    end
    disp(rot_period)

    rot_plot(j) = rot_period;

    j = j+1;
end


%% Plot

[X, Y] = meshgrid(rot_plot, h_plot);
figure
hold on
contour3(X, Y, ratio_v, [30 30], 'r', 'LineWidth', 2);
surf(X, Y, ratio_v)
shading interp
xlabel('Rotation period [h]')
ylabel('Altitude [m]')
zlabel('SR')
colorbar
% clim([0 100])
caxis([0 100])
legend('30% smear')
ylim([1000 15000])