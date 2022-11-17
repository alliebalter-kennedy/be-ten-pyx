function out = LABCO_objective_steadyerosion(X, data, constants, p, dFlag)

% This is the objective function for fitting exercise 2, which is assuming
% steady erosion (rate floats) and a fixed exposure time. 
% out = LABCO_objective_zeroerosion(X,d,p,dFlag)
% 
% where X has the following params:
%   X(1) erosion rate (cm Myr-1)
%   X(2) He-3/Be-10 ratio (non-dimensional)
%   X(3) Lsp (g/cm2) spallogenic e-folding length
%   X(4) k_neg_10 (fraction) negative muon capture cross section
%   X(5) k_neg_3 (fraction) negative muon capture cross section
%   X(6) depthOffset (cm) missing depth below 18 cm
%   X(7) nonCosmoHe3 (atoms g-1) amount of non-cosmogenic He-3
%
% samples is data structure with core data
% p is additional struct with things like surface production rate, etc. 
% dflag is 0 for objective fctn; 1 for diagnostic struct
%
% Modified from script by Greg Balco - Berkeley Geochronology Center - 2018
% Allie Balter-Kennedy - Lamont-Doherty Earth Observatory - 2022
% Not licensed for reuse or distribution

%% Get X in correct dimensions

erosionRate = X(1) .* constants.rho .* 1e-6; % g/cm2/yr
P10sp = X(2).* p.SFsp; % No conversion
Lsp = X(3); % No conversion
k_neg_10 = X(4); % No conversion
k_neg_3 = X(5); % No conversion
depthOffset = X(6); % cm
nonCosmoHe3 = X(7); % Atoms g-1

%% Define production

he.zcm = data.he.depth_cm;
he.zcm(he.zcm>=18) = he.zcm(he.zcm>=18) + depthOffset;
he.zgcm2 = he.zcm .* constants.rho; % depth, g cm-2

be.zcm = data.be.depth_cm;
be.zcm(be.zcm>=18) = be.zcm(be.zcm>=18) + depthOffset;
be.zgcm2 = be.zcm .* constants.rho; % depth, g cm-2

% Spallation

result.N3sp = exp(-he.zgcm2./Lsp).*(p.P3sp.*Lsp./erosionRate).*(1-exp(-p.texp.*(erosionRate./Lsp))); % spallogenic He-3 produced at sample depths
result.N10sp = ((P10sp .* exp(-be.zgcm2./Lsp))./(constants.l10 + erosionRate./Lsp)).*(1-exp(-(constants.l10 + erosionRate./Lsp).*p.texp)); % spallogenic Be-10 produced at sample depths

% muons

% helium-3
for a = 1:length(data.he.N3_normalized)
    result.N3fast(a) = integral(@(t) interp1(log(p.predz),p.m3stub_fast,log(he.zgcm2(a)+(erosionRate.*t))),0,p.texp,'reltol',1e-3,'abstol',1e-1); % He-3 produced by fast muons at sample depths
    result.N3neg(a) = integral(@(t) k_neg_3.*interp1(log(p.predz),p.m3stub_neg,log(he.zgcm2(a)+(erosionRate.*t))),0,p.texp,'reltol',1e-3,'abstol',1e-1); % He-3 produced by negative muon capture at sample depths
end

% beryllium-10
for a = 1:length(data.be.N10)
    result.N10fast(a) = integral(@(t) exp(-constants.l10.*t).*interp1(log(p.predz),p.m10stub_fast,log(be.zgcm2(a)+erosionRate.*t)),0,p.texp,'reltol',1e-3,'abstol',1e-1); % Be-10 produced by fast muons at sample depths
    result.N10neg(a) = integral(@(t) k_neg_10.* exp(-constants.l10.*t).*interp1(log(p.predz),p.m10stub_neg,log(be.zgcm2(a)+erosionRate.*t)),0,p.texp,'reltol',1e-3,'abstol',1e-1); % Be-10 produced by negative muon capture at sample depths
end


% Total
result.N10tot = result.N10sp + result.N10neg' + result.N10fast';
result.N3tot = result.N3sp + result.N3neg' + result.N3fast' + nonCosmoHe3;

% Miss

result.miss10 = (result.N10tot - data.be.N10)./(data.be.dN10); 
result.miss3 = (result.N3tot - data.he.N3_normalized)./(sqrt(data.he.dN3.^2 + 1.1e6.^2));
result.x2 = sum(result.miss10.^2) + sum(result.miss3.^2);

result.depth3 = he.zcm;
result.depth10 = be.zcm;

if dFlag == 0
    out = result.x2;
elseif dFlag == 1
    out = result;
end