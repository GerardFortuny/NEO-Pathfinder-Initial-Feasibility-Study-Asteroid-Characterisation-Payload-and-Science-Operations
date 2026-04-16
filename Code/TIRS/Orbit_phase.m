clear all
close all

%% General data
m = 7.3E+10;          % SU49 asteroid mass [kg]
G = 6.6743E-11;       % Gravitational constant [m^3 kg^-1 s^-2]
mu = G*m;

h = 700;

% OTES
FOV = 6.35e-3; % [rad]
t_int = 1.8; % OTES integration time [s]

%% Minimum time needed to scan the asteroid

j = 1;
for R_a = 200:10:400

    A_T = 4*pi*R_a^2; %Toal area

    D = tan(FOV)*h;
    A_c = pi*(D/2)^2; %Area covered by the instrument

    punts = A_T/A_c;

    temps(j) = (t_int*punts)/3600;

    j = j+1;
end

R_plot = 200:10:400;

%% Resolution
% OTES:
IFOV = FOV;

j = 1;
for R_a = 200:10:400
    lambd = asin((R_a+h)*sin(IFOV/2)/R_a) - (IFOV/2);

    Res(j) = lambd*2*R_a;

    j = j+1;
end

%% Smear
% We are doing polar orbits, so the maximum will be in the equator.

j = 1;
for R_a = 200:10:400
    i = 1;
    for rot_period = 1:0.25:10 % Asteroid rotation period [h]

        v_rot = R_a*2*pi/(rot_period*3600);

        v_sup = sqrt(mu/(R_a+h)^3)*R_a;

        v_T = sqrt(v_rot^2 + v_sup^2);

        smear(j,i) = 100*v_T*t_int/Res(j);

        rot_plot(i) = rot_period;
        i = i+1;


    end

    R_plot(j) = R_a;
    j = j+1;

end


%% Plot

figure
plot(R_plot,Res,'LineWidth',1.5)
xlabel('Asteroid radius [m]')
ylabel('Resolution [m]')
grid on

figure
plot(R_plot,temps,'LineWidth',1.5)
xlabel('Asteroid radius [m]')
ylabel('Time needed [h]')
grid on

[X, Y] = meshgrid(rot_plot, R_plot);
figure
hold on
contour3(X, Y, smear, [10 10], 'r', 'LineWidth', 2);
pcolor(X, Y, smear)
shading interp
xlabel('Rotation period [h]')
ylabel('Asteroid radius [m]')
zlabel('SR')
colorbar
legend('10% smear')





