% This script performs the fitting exercise for the steady erosion case. 

% Allie Balter-Kennedy - Lamont-Doherty Earth Observatory - March 2022
% Utilizes scripts written by Greg Balco, with permission

clear all
close all

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

%% Generate muon fluxes and stopping rates at all sample depths

p.mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
p.mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
% p.mc10.sigma0 = 0.280e-30.*0.5; % ubarns; Balco 2017
% p.mc10.sigma0 = 0.280e-30.*1.5; % ubarns; Balco 2017
p.mc10.k_neg = 1; % Dummy

% p.mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
p.mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
p.mc3.sigma0 = 5.70e-30; % ubarns; from Balco fit to Larsen data
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
k_neg_10_0 = 0.00191 .* 0.704 .*0.1828; % fstar * fC * fD; fC and fD from Heisinger 2002b for O in SiO2. Will fit whole k_neg for pyroxene. 
k_neg_3_0 = 0.0045; % guess from Nesterenok
depthOffset0 = 0; 
nonCosmoHe3_0 = 3.3e6; % from Balco blog post

x0 = [erosionRate0 P10sp0 Lsp0 k_neg_10_0 k_neg_3_0 depthOffset0 nonCosmoHe3_0];
xmin = [0 0 0 0 0 -20 nonCosmoHe3_0];
xmax = [Inf Inf Inf Inf Inf 20 nonCosmoHe3_0];

% fit age only (assume 3 cm missing below 18 cm)
% xmin = [0 P10sp0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
% xmax = [Inf P10sp0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
%% optimization settings

opts = optimset('fmincon');
opts = optimset(opts,'tolfun',1e-7, 'display', 'none');

%% MC simulation

ni = 1000; % number of iterations

for i = 1:ni
% for each iteration, choose values for he, be and P3
    this_p = p;
    this_p.P3sp = (randn(1,1).*(0.108.*constants.P3_SLHL) + constants.P3_SLHL).*p.SFsp; %10.8% uncertainty on P3
    this_data = data;
    this_data.be.N10 = randn(size(this_data.be.N10)) .*this_data.be.dN10 + this_data.be.N10;
    this_data.he.N3_normalized = randn(size(this_data.he.N3_normalized)) .*this_data.he.dN3 + this_data.he.N3_normalized;
    [optx,fval] = fmincon(@(x) LABCO_objective_steadyerosion(x,this_data,constants, this_p,0),x0,[],[],[],[],xmin,xmax,[],opts);
    P10sp(i) = optx(2);
end

%% save 

save('LABCO_MC_steadyerosion.mat', "P10sp")

