clear all
close all

%% Call labco data

[loc he be] = create_LABCO_data();

%% Set up lab flags

ldeo = find(strcmp(he.prep_lab, 'LDEO'));
bgc_bgc = find(strcmp(he.prep_lab, 'BGC') & strcmp(he.analysis_lab, 'BGC'));
bgc_crpg = find(strcmp(he.prep_lab, 'BGC') & strcmp(he.analysis_lab, 'CRPG'));
crpg_crpg = find(strcmp(he.prep_lab, 'CRPG') & strcmp(he.analysis_lab, 'CRPG'));
crpg_bgc = find(strcmp(he.prep_lab, 'CRPG') & strcmp(he.analysis_lab, 'BGC'));

%% Plot He-3 data

figure(1)

% subplot(1, 2, 1)
% hold on
% errorbar(he.N3(ldeo), he.depth_cm(ldeo), he.dN3(ldeo), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
% errorbar(he.N3(bgc_bgc), he.depth_cm(bgc_bgc), he.dN3(bgc_bgc), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
% errorbar(he.N3(bgc_crpg), he.depth_cm(bgc_crpg), he.dN3(bgc_crpg), 'horizontal', 'kv', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
% errorbar(he.N3(crpg_crpg), he.depth_cm(crpg_crpg), he.dN3(crpg_crpg), 'horizontal', 'kv', 'MarkerFaceColor', 'w', 'MarkerSize', 8)
% errorbar(he.N3(crpg_bgc), he.depth_cm(crpg_bgc), he.dN3(crpg_bgc), 'horizontal', 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8)
% set(gca, 'xscale', 'log', 'fontsize', 16, 'ydir', 'reverse', 'xlim', [1e7 1e9], 'ylim', [0 180])
% xlabel('Measured [^{3}He] (atoms g^{-1})')
% ylabel('Depth (cm)')
% % legend('LDEO prep/LDEO meas.', 'BGC prep/BGC meas.', 'BGC prep/CRPG meas.', 'CRPG prep/CRPG meas.', 'CRPG prep/BGC meas.', 'Location', 'Southeast')
% grid on
% box on
% 
% subplot(1, 2, 2)
hold on
errorbar(he.N3_normalized(ldeo), he.depth_cm(ldeo), he.dN3(ldeo), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N3_normalized(bgc_bgc), he.depth_cm(bgc_bgc), he.dN3(bgc_bgc), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N3_normalized(bgc_crpg), he.depth_cm(bgc_crpg), he.dN3(bgc_crpg), 'horizontal', 'kv', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N3_normalized(crpg_crpg), he.depth_cm(crpg_crpg), he.dN3(crpg_crpg), 'horizontal', 'kv', 'MarkerFaceColor', 'w', 'MarkerSize', 8)
errorbar(he.N3_normalized(crpg_bgc), he.depth_cm(crpg_bgc), he.dN3(crpg_bgc), 'horizontal', 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8)

% add exponenetial curve

Lsp = 140; % spallation attenuation length
rho = 2.94; % Ferrar dolerite density
N3_surf = mean(he.N3_normalized(he.depth_cm == round(0.5, 1))) .* exp((0.5.*rho)/Lsp);

yy = [0:180]; % cm

N3_sp = N3_surf .* exp(-(yy.*rho)./Lsp);

hold on
plot(N3_sp, yy, 'k')

set(gca, 'xscale', 'log', 'fontsize', 16, 'ydir', 'reverse', 'xlim', [1e7 1e9], 'ylim', [0 180])
xlabel('Normalized [^{3}He] (atoms g^{-1})')
ylabel('Depth (cm)')
legend('LDEO prep/LDEO meas.', 'BGC prep/BGC meas.', 'BGC prep/CRPG meas.', 'CRPG prep/CRPG meas.', 'CRPG prep/BGC meas.', 'Spallation curve', 'Location', 'Southeast')
grid on
box on

%% Plot Be-10 data

figure(2)


hold on
errorbar(be.N10(be.outlier == 0), be.depth_cm(be.outlier == 0), be.dN10(be.outlier == 0), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(be.N10(be.outlier == 1), be.depth_cm(be.outlier == 1), be.dN10(be.outlier == 1), 'horizontal', 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8)


% add exponenetial curve

N10_surf = mean(be.N10(be.depth_cm == round(24.5, 1))) .* exp((24.5.*rho)/Lsp);
N10_sp = N10_surf .* exp(-(yy.*rho)./Lsp);

hold on 
plot(N10_sp, yy, 'k')

set(gca, 'xscale', 'log', 'fontsize', 16, 'ydir', 'reverse', 'xlim', [1e5 1e7], 'ylim', [0 180])
xlabel('Measured [^{10}Be] (atoms g^{-1})')
ylabel('Depth (cm)')
legend('Measured [^{10}Be]', 'Outlier', 'Spallation curve', 'Location', 'Southeast')
grid on
box on

%% Plot He-4 data

figure(3)

hold on
errorbar(he.N4(ldeo), he.depth_cm(ldeo), he.dN4(ldeo), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N4(bgc_bgc), he.depth_cm(bgc_bgc), he.dN4(bgc_bgc), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N4(bgc_crpg), he.depth_cm(bgc_crpg), he.dN4(bgc_crpg), 'horizontal', 'kv', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(he.N4(crpg_crpg), he.depth_cm(crpg_crpg), he.dN4(crpg_crpg), 'horizontal', 'kv', 'MarkerFaceColor', 'w', 'MarkerSize', 8)
errorbar(he.N4(crpg_bgc), he.depth_cm(crpg_bgc), he.dN4(crpg_bgc), 'horizontal', 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8)
set(gca, 'xscale', 'log', 'fontsize', 16, 'ydir', 'reverse', 'xlim', [0.1e14 3e14], 'ylim', [0 180])
xlabel('Measured [^{4}He] (atoms g^{-1})')
ylabel('Depth (cm)')
legend('LDEO prep/LDEO meas.', 'BGC prep/BGC meas.', 'BGC prep/CRPG meas.', 'CRPG prep/CRPG meas.', 'CRPG prep/BGC meas.', 'Location', 'Southwest')
grid on
box on
