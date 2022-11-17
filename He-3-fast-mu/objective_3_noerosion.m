function out = objective_3_noerosion(X,d,dFlag)

% This is the objective function for muon production model fitting to He-3
% data. 

if nargin < 3; dFlag = 0; end

% This one is stationary, no erosion. 

% Unpack

sigma0 = X(1).*1e-30; % cross-section
alpha = X(2); % alpha
N3nuc = X(3).*1e6; % nucleogenic He-3, atoms/g
pxil = X(4); % conversion between px/ilmenite

t = 16e6; % age

% Compute production rates

P3 = P_mu_3_free_alpha(d.z,929,alpha,sigma0);
%P3 = P_mu_3_free_alpha(d.z,1013.25,alpha,sigma0);


N3p = N3nuc + P3.*t;

c = d.pyroxene.*1 + d.ilmenite./pxil;

N3obs = d.N3raw.*c;

miss = (N3p' - N3obs)./d.dN3raw;

if dFlag == 1
    out.N3p = N3p';
    out.N3obs = N3obs;
    out.miss = miss;
    out.P3 = P3;
else
    out = sqrt(sum(miss.^2));
end



