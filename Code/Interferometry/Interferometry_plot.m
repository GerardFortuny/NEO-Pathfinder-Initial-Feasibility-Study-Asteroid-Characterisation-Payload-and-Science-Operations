clear all
close all

%% Data
x = [100, 1000, 10000];
y = [25, 50, 100, 150, 300, 1000];

Z = [
    5.34e-01    2.942026582    5.024809573;
    0.204778124 1.300465192    2.214005725;
    0.080610173 0.448938562    0.764409627;
    0.035706064 0.199968365    0.343653896;
    0.005939697 0.033700739    0.057384435;
    0.000182434 0.001084493    0.002012535
];


%% Plot

[X, Y] = meshgrid(x, y);

figure
pcolor(X, Y, Z);
shading interp;
hold on;
[C, h] = contour(X, Y, Z, [1 1], 'r', 'LineWidth', 3);
legend(h, 'SNR = 1')
set(gca, 'XScale', 'log', 'YScale', 'log');
colorbar;
xlabel('Eje X (log)');
ylabel('Eje Y (log)');


