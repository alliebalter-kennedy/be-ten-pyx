% this script uses the zero erosion and steady erosion assumptions to solve
% for the maximum and minimum be-10 p rates in pyroxene.

clear all
close all


%% Add relevant paths and load data

addpath('data/')

filename = 'SULA_data.txt'; % file where sample, He-3, and Be-10 data are stored.
data = readtable(filename); % load data

%% Unpack data into usable form

% put all sample information into structural array for use below.

% put Sample IDs in array
samples.ID = table2cell(data(:, 'Sample_ID')); % Sample IDs

% put sample data in array
samples.thickness = table2array(data(:, 'thickness')); % Sample thickness in g cm-2

% put He-3 data in array
samples.N3_meas = table2array(data(:, 'N3_LDEO')); % He-3 concentrations
samples.dN3 = table2array(data(:, 'dN3_LDEO')); % and uncertainty, measured at LDEO, [atoms g^-1]

                                                      
% CRONUS-P (CPX-2) He-3 concentration measured at LDEO; [atoms ^g-1]
samples.cronusPMeasured = zeros(5, 1);
samples.cronusPMeasured(:) = samples.N3_meas(strcmp(samples.ID, 'CPX-2'));

% normalize He-3 concentrations

samples.N3 = samples.N3_meas .* (5.02E+09./samples.cronusPMeasured); 

% put Be-10 data in array
samples.N10 = table2array(data(:, 'N10'));      % Be-10 concentrations
samples.dN10 = table2array(data(:, 'dN10'));    % and uncertainty; prepared
                                                % at LDEO, measured at
                                                % LLNL-CAMS; [atoms g^-1]

% put location data in array
samples.lat = table2array(data(:, 'lat'));      % Latitude; [DD]
samples.long = table2array(data(:, 'long'));    % Longitude; [DD]
samples.elv = table2array(data(:, 'elv'));      % Elevation; [m]
samples.atm = table2cell(data(:, 'atm'));       % Amosphere model; 'ant' is
                                                % Antarctic; 'std' is
                                                % standard
p.N3noncosmo = 3.3e6; % atoms/g
p.dN3noncosmo = 1.1e6;
%% Constants 

p.constants = bedrock_constants();

p.constants.dP3p_St = 0.108 .* p.constants.P3p_St;

% other constants

p.constants.rho = 2.94;
p.constants.Lsp = 140;

% change sample thickness to g cm-2

samples.thickness = samples.thickness .* p.constants.rho;

% Atmospheric pressure at site
p.h = antatm(samples.elv); % site pressure

% Define p rate info
p.sf_St = stone2000(samples.lat,p.h,1); % scaling factor

p.sf_thick = (p.constants.Lsp./samples.thickness) .* (1-exp(-(samples.thickness./p.constants.Lsp)));

p.sf_tot = p.sf_St' .* p.sf_thick;

%% muon set up

mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
mc10.k_neg = 0.000174; % from SULA steady erosion

% mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
mc3.sigma0 = 6.01e-30; % ubarns; from Balco fit to Larsen data
mc3.k_neg = 0.00364; % from SULA steady erosion

for a = 1:length(p.h)   
   P3mu(a) = P_mu_total_alpha1((samples.thickness(a)./2),p.h(a),mc3);
end

p.texp = 14.5.*1e6; 

%% solve cases

% zero erosion -- solve for t and P10(0) for each sample

% first, get t using spallation and muons

tt = (samples.N3 - p.N3noncosmo)./((p.constants.P3p_St.*p.sf_tot) + P3mu');

for a = 1:length(p.h)
    P10mu(a, 1) = P_mu_total_alpha1((samples.thickness(a)./2),p.h(a),mc10);
end

N10_mu_tt = (P10mu ./ p.constants.l10) .* (1-exp(-p.constants.l10.*tt));

P10_tt = ((samples.N10-N10_mu_tt) .* p.constants.l10) ./ (p.sf_tot .* (1-exp(-p.constants.l10.*tt)));

%% steady erosion -- solve for ee and P10(0)


for a = 1:length(p.h)
    this_p.constants = p.constants;
    this_p.h = p.h(a);
    this_thick = samples.thickness(a);
    this_p.P3sp = p.constants.P3p_St.*p.sf_tot(a);
    this_N3 = samples.N3(a);
    this_dN3 = samples.dN3(a);
    this_p.texp = p.texp;
    this_p.N3noncosmo = p.N3noncosmo;
    this_p.dN3noncosmo = p.dN3noncosmo;
    
    x_0 = 5 .* 1e-6 .* p.constants.rho; % 5 cm Myr-1,converted to g/cm2 per yr initial guess
    ee(a) = fminsearch(@(x) get_erosion(x, this_N3, this_dN3, this_thick, this_p, mc3), x_0);
end

opt_ee_cm = ee ./ (1e-6 .* p.constants.rho);
%%
for a = 1:length(p.h)
    N10_mu_ee(a, 1) = integral(@(t) exp(-p.constants.l10.*t).* P_mu_total_alpha1((samples.thickness(a)./2) + (ee(a) .* t), p.h(a), mc10), 0, p.texp, "RelTol", 1e-3);
end

P10_ee = (samples.N10 - N10_mu_ee).* (p.constants.l10 + (ee'./p.constants.Lsp))./p.sf_tot;
    
%% uncertainties

dP10_tt = sqrt((samples.dN3./samples.N3_meas).^2 + (samples.dN10./samples.N10).^2 + (p.constants.dP3p_St./p.constants.P3p_St).^2) .* P10_tt;
dP10_ee = sqrt((samples.dN3./samples.N3_meas).^2 + (samples.dN10./samples.N10).^2 + (p.constants.dP3p_St./p.constants.P3p_St).^2) .* P10_ee;

%% save

save('SULA_minmax.mat', "P10_tt", "P10_ee", 'dP10_ee', 'dP10_tt', 'samples');
    
%% Objective function for fitting erosion rate


function out = get_erosion(X, measured_N3, dN3, thick, p, mc3)

erosionRate = X(1); 

pred.N3mu = integral(@(t) P_mu_total_alpha1((thick./2) + (erosionRate.* t),p.h,mc3), 0, p.texp, 'reltol',1e-3);

pred.N3sp = (p.P3sp.*p.constants.Lsp./erosionRate).*(1-exp(-p.texp.*(erosionRate./p.constants.Lsp)));

pred.N3 = pred.N3mu + pred.N3sp + p.N3noncosmo;

x2 = ((pred.N3 - measured_N3)./(sqrt(dN3.^2 + p.dN3noncosmo.^2))).^2;

out = x2;
end