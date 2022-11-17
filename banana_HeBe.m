% This wrapper script creates the two-nuclide diagram (Figure 8). The script
% banana_background.m is required to generate the background image.

clear all
close all

t = [0:1e3:50e6];
erosion = [0:1e-6:1000e-5];

l10 = 4.9975E-07;
rho = 2.94;
Lsp = 140;

constants.cronusPAccepted = 5.02E+09;   % Accepted CRONUS-P concentration; 
                                        % Blard et al. 2015
constants.N3_nonCosmogenic = 3.3e6;       % Concentration of non-cosmogenic 
                                        % He-3 in Ferrar Dolerite; [atoms
                                        % g^-1]; Balter-Kennedy et al.,
                                        % 2020 and references therein
%% Load LABCO data
% Location data

loc.latitude = -77.54976;       % [DD]
loc.longitude = 160.9578;       % [DD]
loc.elevation = 990.2;          % [m]

loc.atm = 'ant';                % Atmosphere model flag for v3
loc.shielding = 1;              % shielding
loc.yearCollected = 2010;       % year core was collected


%% Add relevant paths

addpath('data/')

%% Load data

% Note that numbered surface samples in Balter-Kennedy, submitted Nov. 2022
% are the same as the samples listed here with the prefix SULA. 

labco_filename = 'LABCO_data.txt';

labco_data = readtable(labco_filename);             % load data

sula_filename = 'SULA_data.txt';
sula_data = readtable(sula_filename);  

%% Unpack data into usable form

% put all sample information into structural array for use below.

% put Sample IDs in array
sula.ID = table2cell(sula_data(:, 'Sample_ID')); % Sample IDs

% put sample data in array
sula.thickness = table2array(sula_data(:, 'thickness')); % Sample thickness in g cm-2

% put He-3 data in array
sula.N3_meas = table2array(sula_data(:, 'N3_LDEO')); % He-3 concentrations
sula.dN3 = table2array(sula_data(:, 'dN3_LDEO')); % and uncertainty, measured at LDEO, [atoms g^-1]

                                                      
% CRONUS-P (CPX-2) He-3 concentration measured at LDEO; [atoms ^g-1]
sula.cronusPMeasured = zeros(5, 1);
sula.cronusPMeasured(:) = sula.N3_meas(strcmp(sula.ID, 'CPX-2'));


% put Be-10 data in array
sula.N10 = table2array(sula_data(:, 'N10'));      % Be-10 concentrations
sula.dN10 = table2array(sula_data(:, 'dN10'));    % and uncertainty; prepared
                                                % at LDEO, measured at
                                                % LLNL-CAMS; [atoms g^-1]

% put location data in array
sula.lat = table2array(sula_data(:, 'lat'));      % Latitude; [DD]
sula.long = table2array(sula_data(:, 'long'));    % Longitude; [DD]
sula.elv = table2array(sula_data(:, 'elv'));      % Elevation; [m]
sula.atm = table2cell(sula_data(:, 'atm'));       % Amosphere model; 'ant' is
                                                % Antarctic; 'std' is
                                                % standard

%% add SULA-444 - not really efficient
sula.ID(end+1) = {'SULA-444'};
sula.thickness(end+1) = 1.5;
sula.N3_meas(end+1) = 5.20E+08;
sula.dN3(end+1) = 2.00E+07;
sula.cronusPMeasured(end+1) = 5.21E+09;
sula.N10(end+1) = 10204465;
sula.dN10(end+1) = 245515;
sula.lat(end+1) = -77.52;
sula.long(end+1) = 160.91;
sula.elv(end+1) = 1145;
sula.atm(end+1) = {'ant'};


%% add LABCO, same, not efficient 
sula.ID = {'318', '439', '446s', '464', 'NXP 93*52', '444'}; % transform sample names so same as text

sula.ID(end+1) = {'LABCO'};
sula.thickness(end+1) = 1;
sula.N3_meas(end+1) = 4.54E+08;
sula.dN3(end+1) = 9.23E+06;
sula.cronusPMeasured(end+1) = 4.78E+09;
sula.N10(end+1) = 8.29e6;
sula.dN10(end+1) = 8.29e6.*0.03;
sula.lat(end+1) = -77.54976;
sula.long(end+1) = 160.9578;
sula.elv(end+1) = 990.2;
sula.atm(end+1) = {'ant'};
%% normalize He-3 concentrations

sula.N3 = sula.N3_meas .* (5.02E+09./sula.cronusPMeasured); 
%% Production rate info

p.constants = bedrock_constants();

p.constants.P10p_St = 3.6; % from summary stats exercise in section 5.1.3
p.constants.P3p_St = 124.03; % from Borchers et al., 2016
HeBeRatio = p.constants.P3p_St./p.constants.P10p_St; 

% Atmospheric pressure at site
p.siteP = antatm(sula.elv); % site air pressure

% Define St scaling factor
p.SFsp = stone2000(sula.lat,p.siteP,1); % scaling factor

p.sf_thick = (Lsp./sula.thickness) .* (1-exp(-(sula.thickness./Lsp)));

% Get production rates

p.P3sp = p.constants.P3p_St .* p.SFsp' .* p.sf_thick;
p.P10sp = p.constants.P10p_St .* p.SFsp' .* p.sf_thick;

%% Normalize Be and He concentrations to PR

banana.HeNorm = (sula.N3 - constants.N3_nonCosmogenic) ./ p.P3sp;
banana.BeNorm = sula.N10 ./ p.P10sp;

banana.dHeNorm = sula.dN3./ p.P3sp;
banana.dBeNorm = sula.dN10./ p.P10sp;
%% plot

[bananaPlot ax1] = banana_background_simple(HeBeRatio, Lsp, l10, rho, erosion, t);

%%
axes(ax1)
hold on

for a = 1:length(sula.dN3)
    [x y] = ellipse(banana.HeNorm(a), banana.dHeNorm(a), banana.BeNorm(a), banana.dBeNorm(a), 0);

    fill(x, y, [0.9 0.9 0.9])
    plot(x, y, 'k')
    if (mean(y) >= 0.3 & mean(y) <= 0.4) & strcmp(sula.ID(a), 'LABCO') == 0
        x_temp = max(x(find((abs(y - (mean(y)+(0.068./2))) == min(abs(y - (mean(y)+(0.068./2))))))));
        text(x_temp - 30000, mean(y), sula.ID(a), 'BackgroundColor', [1 1 1], 'Margin', 0.01, 'HorizontalAlignment','right', 'FontSize', 12)
    elseif (mean(y) > 0.4 | mean(y) < 0.3) & strcmp(sula.ID(a), 'LABCO') == 0
        x_temp = max(x(find((abs(y - (mean(y)-(0.068./2))) == min(abs(y - (mean(y)-(0.068./2))))))));
        text(x_temp + 20000, mean(y), sula.ID(a), 'BackgroundColor', [1 1 1], 'Margin', 0.01, 'HorizontalAlignment', 'left', 'FontSize', 12)
    elseif strcmp(sula.ID(a), 'LABCO') == 1
        x_temp = max(x(find((abs(y - (mean(y)+(0.068./2))) == min(abs(y - (mean(y)+(0.068./2))))))));
        text(x_temp - 30000, mean(y), sula.ID(a), 'BackgroundColor', [1 1 1], 'Margin', 0.01, 'HorizontalAlignment','right', 'FontSize', 12)

    end
end
plot(banana.HeNorm, banana.BeNorm./banana.HeNorm, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)


% axes(ax2)
% hold on
% for a = 1:length(sula.dN3)
%     [x y] = ellipse(banana.HeNorm(a), banana.dHeNorm(a), banana.BeNorm(a), banana.dBeNorm(a), 0);
%     fill(x, y, [0.9 0.9 0.9])
%     plot(x, y, 'k')
% end
% 
% plot(banana.HeNorm, banana.BeNorm./banana.HeNorm, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 2)
