% Relative error estimate for Be-10 in pyroxene given existing measurements. This is really greg balcos sript

%% Clean up

clear all; close all;

%% Load Be measurement compilations and determine concentration-error 
% relationship

summary_Be_errors = [5.06E+07	1.21E+06
1.55E+07	3.72E+05
2.84E+07	6.80E+05
2.75E+07	6.57E+05
6.46E+07	1.54E+06
1.25E+07	2.96E+05
1.06E+07	2.53E+05
1.07E+07	3.38E+05
1.04E+07	2.51E+05
1.02E+07	2.46E+05
5.21E+06	1.45E+05
4.39E+06	1.17E+05
3.28E+06	1.06E+05
2.60E+06	8.11E+04
1.57E+06	5.90E+04
9.01E+05	4.59E+04
1.11E+06	4.32E+04
7.78E+05	4.22E+04
5.25E+05	3.54E+04
4.32E+05	3.86E+04];

pyx_wts = [0.0996
0.1113
0.1019
0.1048
0.1150
0.1537
0.2228
0.1971
0.1949
0.1571
0.1486
0.1500
0.1510
0.1472
0.1246
0.1468
0.1801
0.1225
0.1412
0.1109];

N10 = summary_Be_errors(:,1);
delN10 = summary_Be_errors(:,2);
rd10 = delN10./N10;

% Curve fits

% Log-linear

p10 = polyfit(log(N10),log(delN10),1);

plotN10 = logspace(3,9,20);
plotdelN10 = exp(p10(1).*log(plotN10) + p10(2));
plotrd10 = plotdelN10./plotN10;


p10c = polyfit(N10.*pyx_wts, delN10.*pyx_wts, 1);
del10c = p10c(1).*(N10.*pyx_wts) + p10c(2);

figure
plot(N10.*pyx_wts, delN10.*pyx_wts, 'r.')
hold on
plot(N10.*pyx_wts, del10c, 'r');

plotlin = [1e3:1e3:1e8];
del10cplot = p10c(1).*plotlin + p10c(2); 
rd10c = del10cplot./plotlin;

figure
plot(N10.*pyx_wts, rd10, 'r.')
hold on
plot(plotlin, rd10c, 'r')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
set(gca, 'xlim', [1e3 1e8]);
grid on
set(gca, 'fontsize', 14)
xlabel('^{10}Be atoms', 'Fontsize', 16)
ylabel('Relative uncertainty', 'Fontsize', 16)

