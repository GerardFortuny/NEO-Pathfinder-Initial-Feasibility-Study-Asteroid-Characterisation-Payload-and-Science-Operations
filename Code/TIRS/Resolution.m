clear all
close all

%% General data
FOV = 6.35e-3; % [rad]
R_a = 300; % Asteroid radius [m]
h = 4000;

%% Calculations

lambd_max = acos(R_a/(R_a+h));

D = [];
q = 1;
for d = 0:0.01:cos(pi/2 - lambd_max)*R_a

    r = sqrt(R_a^2 - d^2);
    lambd_s = atan2(d,r);
    s_s = sqrt(R_a^2 + (R_a + h)^2 -2*R_a*(R_a+h)*cos(lambd_s));
    eta_s = asin(R_a*sin(lambd_s)/s_s);
    del_lambd = lambd_s - asin((R_a+h)*sin(eta_s - FOV)/R_a) + (eta_s - FOV);

    lambd_c = asin((R_a+h)*sin(eta_s - FOV/2)/R_a) - (eta_s - FOV/2);

    if lambd_c > 0
        D(q,1) = del_lambd*R_a;
        D(q,2) = lambd_s;

        Cov(q) = sin(lambd_s)*100;
        q = q+1;

    end
end

%% Plot

figure
plot(Cov, D(:,1),LineWidth=1.5)
yline(31,'r',LineWidth=1.5)
xlabel('Covered part of the path [%]')
ylabel('Resolution [m]')
legend('Achieved resolution', 'Required resolution')
grid on
ylim([20 130])

