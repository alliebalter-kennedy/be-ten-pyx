% This script performs the fitting exercise for the steady erosion case. 

% Allie Balter-Kennedy - Lamont-Doherty Earth Observatory - March 2022
% Utilizes scripts written by Greg Balco, with permission

clear all
close all

addpath('data/')

%% Define constants

%% constants

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator

constants.l10 = 4.9975E-07;             % Be-10 decay rate; [yr^-1];
                                        % Nishizumii 2007? Check reference.

constants.rho = 2.94;                   % Density of Ferrar Dolerite; 
                                        % [g cm^-3]

constants.Lsp = 140;                    % spallation attenuation length for 
                                        % sample thickness corrections.
                                        % Should be appropriate for
                                        % Antarctica.

p.texp = 14.5e6;


%% Load data 

addpath('data/')
[loc data.he data.be] = create_LABCO_data('ignore_outliers');

%% surface production

p.h = antatm(loc.elv); % site pressure
p.SFsp = stone2000(loc.lat,p.h,1); % scaling factor

p.P3sp = constants.P3_SLHL.*p.SFsp;

%% Generate muon fluxes and stopping rates at all sample depths

p.mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
p.mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
% p.mc10.sigma0 = 0.280e-30.*0.5; % ubarns; Balco 2017
% p.mc10.sigma0 = 0.280e-30.*1.5; % ubarns; Balco 2017
p.mc10.k_neg = 1; % Dummy

% p.mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
p.mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
p.mc3.sigma0 = 6.01e-30; % ubarns; from Balco fit to Larsen data
% p.mc3.sigma0 = 5.70e-30.*0.5; % ubarns; from Balco fit to Larsen data
% p.mc3.sigma0 = 5.70e-30.*1.5; % ubarns; from Balco fit to Larsen data
p.mc3.k_neg = 1; % Dummy

p.predz = logspace(0,4,100); %[g cm^-2]

% muon fluxes & stopping rate for He-3 and Be-10
for a = 1:length(p.predz)
    % helium
    this_mu3 = P_mu_total_alpha1(p.predz(a),p.h,p.mc3,'yes');
    p.m3stub_fast(a) = this_mu3.P_fast;
    p.m3stub_neg(a) = this_mu3.P_neg;
end

for a = 1:length(p.predz)
    % beryllium 
    this_mu10 = P_mu_total_alpha1(p.predz(a),p.h,p.mc10,'yes');
    p.m10stub_fast(a) = this_mu10.P_fast;
    p.m10stub_neg(a) = this_mu10.P_neg;
end;


%% Initial guess

% Initial guess 

erosionRate0 = 5; % cm/Myr
P10sp0 = 3.6; % from Eaves et al., 2018
Lsp0 = 140; % g/cm2
% k_neg_10_0 = 0.00191 .* 0.704 .*0.1828; % fstar * fC * fD; fC and fD from Heisinger 2002b for O in SiO2. Will fit whole k_neg for pyroxene. 
k_neg_10_0 = 0.00191 .* 0.520 .*0.1828; % with f_c calulated for 10Be production on O in pyroxene. 
k_neg_3_0 = 0.0045; % guess from Nesterenok
depthOffset0 = 0; 
nonCosmoHe3_0 = 3.3e6; % from Balco blog post

x0 = [erosionRate0 P10sp0 Lsp0 k_neg_10_0 k_neg_3_0 depthOffset0 nonCosmoHe3_0];
xmin = [0 0 0 0 0 -20 nonCosmoHe3_0];
xmax = [Inf Inf Inf Inf Inf 20 nonCosmoHe3_0];

% fit age only (assume 3 cm missing below 18 cm)
% xmin = [0 P10sp0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
% xmax = [Inf P10sp0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
%% Try optimizing

opts = optimset('fmincon');
opts = optimset(opts,'tolfun',1e-7);

[optx,fval] = fmincon(@(x) LABCO_objective_steadyerosion(x,data,constants, p,0),x0,[],[],[],[],xmin,xmax,[],opts);

opt.erosionRate = optx(1);
opt.P10sp = optx(2);
opt.Lsp = optx(3);
opt.k_neg_10 = optx(4);
opt.k_neg_3 = optx(5);
opt.depthOffset = optx(6);
opt.nonCosmoHe3 = optx(7);

% Report results
disp(['Erosion rate = ' sprintf('%0.2f',optx(1)) ' cm Myr^{-1}; P10sp = ' sprintf('%0.2f',optx(2))]);
disp(['Lsp = ' sprintf('%0.2f',optx(3)) ' g cm^{-2}; k_neg_10 = ' sprintf('%0.5f',optx(4)) '; k_neg_3 = ' sprintf('%0.5f',optx(5))]);
disp(['Depth missing = ' sprintf('%0.2f',optx(6)) ' cm;']);

%% Plot results

% Get diagnostics
optd = LABCO_objective_steadyerosion(optx,data,constants, p,1);

pz = [0:5:275].*constants.rho;

mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
% mc10.sigma0 = 0.280e-30.*0.5; % ubarns; Balco 2017
% mc10.sigma0 = 0.280e-30.*1.5; % ubarns; Balco 2017
mc10.k_neg = opt.k_neg_10; % Dummy

% mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
mc3.sigma0 = 6.01e-30; % ubarns; from Balco fit to Larsen data
% mc3.sigma0 = 5.70e-30.*0.5; % ubarns; from Balco fit to Larsen data
% mc3.sigma0 = 5.70e-30.*1.5; % ubarns; from Balco fit to Larsen data
mc3.k_neg = opt.k_neg_3; % Dummy


% calculate predicted depth profiles
    %spallation
        pred.N10sp = (opt.P10sp.*p.SFsp.*exp(-pz./opt.Lsp))./(constants.l10 + ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* (1-exp(-(constants.l10 + ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* p.texp));
        pred.N3sp = (p.P3sp.*exp(-pz./opt.Lsp)./((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* (1-exp(-p.texp .* ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)));

    %muons
        for a = 1:length(pz);
            pred.N10mu(a) = integral(@(t) exp(-constants.l10.*t).*P_mu_total_alpha1((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc10),0,p.texp,'reltol',1e-3);
            pred.N10mu_fast(a) = integral(@(t) exp(-constants.l10.*t).*P_mu_total_alpha1_split((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc10, 'fast'),0,p.texp,'reltol',1e-3);
            pred.N10mu_neg(a) = integral(@(t) exp(-constants.l10.*t).*P_mu_total_alpha1_split((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc10, 'neg'),0,p.texp,'reltol',1e-3);
            pred.N3mu(a) = integral(@(t) P_mu_total_alpha1((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc3),0,p.texp,'reltol',1e-3);
            pred.N3mu_fast(a) = integral(@(t) P_mu_total_alpha1_split((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc3, 'fast'),0,p.texp,'reltol',1e-3);
            pred.N3mu_neg(a) = integral(@(t) P_mu_total_alpha1_split((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc3, 'neg'),0,p.texp,'reltol',1e-3);
        end;

pred.N10 = pred.N10sp + pred.N10mu;


pred.N3tot = pred.N3sp + pred.N3mu + opt.nonCosmoHe3;
pred.N3cosmo = pred.N3sp + pred.N3mu;


%% Plots

% set up
ldeo = find(strcmp(data.he.prep_lab, 'LDEO'));
bgc = find(strcmp(data.he.prep_lab, 'BGC'));
bgc_crpg = find(strcmp(data.he.prep_lab, 'BGC') & strcmp(data.he.analysis_lab, 'CRPG'));

lsidex = 0.1;
lsidew = 0.53;
rsidex = 0.68;
rsidew = 0.28;
ht1 = 0.4;
bot1 = 0.08;
ht2 = 0.4;
bot2 = bot1 + ht1 + 0.1; 

ax1 = axes('pos',[lsidex bot1 lsidew ht1]);
ax11 = axes('pos',[rsidex bot1 rsidew ht1]);

ax2 = axes('pos',[lsidex bot2 lsidew ht2]);
ax21 = axes('pos',[rsidex bot2 rsidew ht2]);

% new sample depths 

data.be.depth_fit = data.be.depth_cm;
data.he.depth_fit = data.he.depth_cm;

data.be.depth_fit(data.be.depth_fit>=18) = data.be.depth_fit(data.be.depth_fit>=18) + opt.depthOffset;
data.he.depth_fit(data.he.depth_fit>=18) = data.he.depth_fit(data.he.depth_fit>=18) + opt.depthOffset;

%% Plot concentrations

% beryllium
figure(1)
axes(ax2)
hold on
plot(pred.N10mu, pz./constants.rho, '-', 'Color', [0.5 0.5 0.5])
plot(pred.N10mu_fast, pz./constants.rho, '-.', 'Color', [0.5 0.5 0.5])
plot(pred.N10mu_neg, pz./constants.rho, '--', 'Color', [0.5 0.5 0.5])
plot(pred.N10sp, pz./constants.rho, 'k--')
plot(pred.N10, pz./constants.rho, 'k')
errorbar(data.be.N10, data.be.depth_fit, data.be.dN10, 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)

set(gca, 'YDir', 'reverse', 'xscale', 'log', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('[^{10}Be]', 'FontSize', 16)
ylabel('Depth (cm)', 'FontSize', 16)
grid on
box on


txt = {['e-rate = ' sprintf('%0.2f',optx(1)) ' cm Myr^{-1}'], ...
    ['P_{10, sp} = ' sprintf('%0.2f',optx(2)) ' at. g^{-1} yr^{-1}'], ['L_{sp} = ' sprintf('%0.2f',optx(3)) ' g cm^{-2}'] ...
    ['f^{*}_{10} = ' sprintf('%0.5f',optx(4)./(0.520 .* 0.1828))], ['f^{*}_{3}f_{C}f_{D} = ' sprintf('%0.5f',optx(5))], ...
    ['Depth offset = ' sprintf('%0.2f',optx(6)) ' cm'], ...
    ['\chi^{2} = ' sprintf('%0.1f', fval)]};

text(1.05e6, 190, txt, 'FontSize', 14)

axes(ax21)
hold on
plot([0 0], [0 pz(end)./constants.rho], 'k-')
errorbar(optd.miss10, data.be.depth_fit, ones(length(data.be.depth_fit), 1), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)

set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('Normalized residual', 'FontSize', 16)
yticklabels([])
grid on
box on

% helium
axes(ax1)
hold on
errorbar(data.he.N3_normalized(ldeo), data.he.depth_fit(ldeo), data.he.dN3(ldeo), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(data.he.N3_normalized(bgc), data.he.depth_fit(bgc), data.he.dN3(bgc), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(data.he.N3_normalized(bgc_crpg), data.he.depth_fit(bgc_crpg), data.he.dN3(bgc_crpg), 'horizontal', 'kv', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
plot(pred.N3mu, pz./constants.rho, '-', 'Color', [0.5 0.5 0.5])
plot(pred.N3mu_fast, pz./constants.rho, '-.', 'Color', [0.5 0.5 0.5])
plot(pred.N3mu_neg, pz./constants.rho, '--', 'Color', [0.5 0.5 0.5])
plot(pred.N3sp, pz./constants.rho, 'k--')
plot(pred.N3cosmo, pz./constants.rho, 'k:')
plot(pred.N3tot, pz./constants.rho, 'k')

% plot(optd.N3tot, zcm, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6)

set(gca, 'YDir', 'reverse', 'xscale', 'log', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('[^{3}He]', 'FontSize', 16)
ylabel('Depth (cm)', 'FontSize', 16)
grid on
box on
legend('LDEO prep/LDEO meas', 'BGC prep/BGC meas', 'BGC prep/CRPG meas', 'Total muons', 'Fast muons', 'Negative muons', 'Spallation', 'Total cosmogenic ^{3}He', 'Total ^{3}He/cosmogenic ^{10}Be', 'Location', 'Southeast')

axes(ax11)
hold on
plot([0 0], [0 pz(end)./constants.rho], 'k-')
errorbar(optd.miss3(ldeo), data.he.depth_fit(ldeo), ones(length(ldeo), 1), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(optd.miss3(bgc), data.he.depth_fit(bgc), ones(length(bgc), 1), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(optd.miss3(bgc_crpg), data.he.depth_fit(bgc_crpg), ones(length(bgc_crpg), 1), 'horizontal', 'kv', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)

set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('Normalized residual', 'FontSize', 16)
yticklabels([])
grid on
box on
