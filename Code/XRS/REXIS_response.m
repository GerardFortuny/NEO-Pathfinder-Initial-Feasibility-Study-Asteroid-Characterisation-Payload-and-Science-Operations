clear all
close all

%% Data
Peaks = [276,339,447,521,567,635,703,846,1040,1245,1484,1735,2300,3303,...
        3685,4004,4500,4916,5401,5885,5936,6387,6478,7048,7453;

        2.8883114,0.4340315,0.1160281,12361.406,1.0203132,0.7721222,...
        2000,1.9191448,43.535076,632.25888,30.022212,329.95207,...
        9.3120083,0.0702551,0.5526232,0.0846011,0.0062748,0.0009606,...
        0.0068857,0.0018753,0.0010541,0.0493573,0.0002871,0.0072805,...
        0.0002173];

scatt_fun = scatter_fun;
Q_fun = Q_values;
internal_noise = internal_noise_fun;

delta_V = 15; % eV
rang = 500:delta_V:6500;

G = 4.24;
F = 0.405;
T_obs = 100*3600; % [s]
T_cal = 100*3600; % [s]
FWHM = 100; % Mean at -60º [eV]


%% Calculations
%Initialisation:
I_bin = zeros(size(rang));

% We obtain C1:
for i = 1:length(rang)
    
    delta_kV_peaks = 1/1000; % Because our data is /keV. We consider the peaks 1 eV wide
    delta_kV = delta_V/1000;

    idx = find(Peaks(1,:) >= rang(i) & Peaks(1,:) < rang(i)+delta_V); % We search for the peaks inside the considered bin
    xmid = rang(i) + delta_V/2;

    I_bin(i) = sum(Peaks(2,idx))*delta_kV_peaks + scatt_fun(xmid)*delta_kV;

    C1(i) = I_bin(i)*G*T_obs*F*Q_fun(xmid); % Signal

    B1(i) = internal_noise(xmid) * T_obs * delta_kV * Q_fun(xmid); % Noise from observation

    A1(i) = internal_noise(xmid) * T_cal * delta_kV * Q_fun(xmid); % Noise from calibration

end

% We apply poisson to obtain C2:
B2 = poissrnd(B1);
C2 = poissrnd(C1) + B2;

A2 = poissrnd(A1);

% Gauss

sigma_eV = FWHM / 2.355; % Calculate the standard deviation of the gausian with the FWHM [https://mathworld.wolfram.com/GaussianFunction.html]

sigma_bins = sigma_eV / delta_V; % Convert the sigma from eV to "bins"

x = -5*sigma_bins : 5*sigma_bins;
Gauss_distr = exp(-(x.^2)/(2*sigma_bins^2)); 
Gauss_distr_norm = Gauss_distr / sum(Gauss_distr); % Normalise the distribution

C3 = conv(C2, Gauss_distr_norm, 'same'); % Convolute, so each bin recives the contribution of the others


% Calculate subtracted signal

C3_R = C3 - A2*T_obs/T_cal; % Noise subtraction


%% Plot
% To avoid negative or very low values for the plot
for i = 1:length(rang)
    if C3(i) <= 1e-2
        C3(i) = 1e-2;
    end
    if C3_R(i) <= 1e-2
        C3_R(i) = 1e-2;
    end
end


figure
stairs(rang, C3_R,'k','LineWidth',1)
set(gca,'YScale','log')
x1 = xline(705,'LineWidth',1,'Color',[1 0 0]);       % Fe-L (vermell)
x2 = xline(1040.98,'LineWidth',1,'Color',[0 0.5 1]); % Na-Ka (blau clar)
x3 = xline(1253.6,'LineWidth',1,'Color',[0 0.8 0]);  % Mg-Ka (verd)
x4 = xline(1486.7,'LineWidth',1,'Color',[0.8 0 0.8]);% Al-Ka (magenta)
x5 = xline(1739.98,'LineWidth',1,'Color',[1 0.6 0]); % Si-Ka (taronja)
x6 = xline(2307.84,'LineWidth',1,'Color',[0 0.9 0.9]);% S-Ka (cian)
x7 = xline(3691.68,'LineWidth',1,'Color',[0.6 0.2 0]);% Ca-Ka (marró)
legend([x1 x2 x3 x4 x5 x6 x7], ...
    {'Fe-L','Na-K\alpha','Mg-K\alpha','Al-K\alpha','Si-K\alpha','S-K\alpha','Ca-K\alpha'})
grid on
xlim([665 6500])
ylim([1e0 inf])
xlabel('Energy (eV)')
ylabel('Counts')


