% this script uses the zero erosion and steady erosion assumptions to solve
% for the maximum and minimum be-10 p rates in pyroxene.

clear all
close all

%% other script below here. 

% This script is supposed to use min/max constraints to get a summary
% distribution. 
% 

rs = load('SULA_minmax_w444.mat');

rc_min =  load('LABCO_MC_zeroerosion.mat');

rc_max = load('LABCO_MC_steadyerosion.mat');

rc.min = mean(rc_min.P10sp);
rc.dmin = std(rc_min.P10sp);

rc.max = mean(rc_max.P10sp);
rc.dmax = std(rc_max.P10sp);


% Define data

mins = [rs.P10_tt; rc.min];
maxs = [rs.P10_ee; rc.max];

% Define uncertainties

dmins = [rs.dP10_tt; rc.dmin];
dmaxs = [rs.dP10_ee; rc.dmax];

% Make a plot 

figure; subplot(2,1,1);
for a = 1:length(mins)
    % Make boxes showing individual uncerts
    xx = [1 -1 -1 1 1].*dmins(a) + mins(a);
    yy = [1 1 -1 -1 1].*0.1 + (length(mins)+1 - a);
    fill(xx,yy,[0.9 0.9 0.9],'edgecolor','none'); hold on;
    xx = [1 -1 -1 1 1].*dmaxs(a) + maxs(a);
    yy = [1 1 -1 -1 1].*0.1 + (length(mins)+1 - a);
    fill(xx,yy,[0.9 0.9 0.9],'edgecolor','none');
    % Dots and connector lines
    xx = [mins(a) maxs(a)];
    yy = [1 1].*(length(mins)+1 - a);
    plot(xx,yy,'k','linewidth',0.5);
    plot(mins(a),(length(mins)+1 - a),'>k','markerfacecolor','k','markersize',8);
    plot(maxs(a),(length(mins)+1 - a),'<k','markerfacecolor','k','markersize',8);
end

ids = {'318', '439', '446s', '464', 'NXP 93*52', '444', 'LABCO'};

set(gca,'xlim',[2 6], 'ylim', [0 length(maxs)+1],'FontSize', 16, 'ytick',[1:7],'yticklabel',flip(ids));
grid on;
box on;
    

%% Now try a Monte Carlo simulation to come up with a permissible
% distribution for the right value. 

ni = 100000;

xgrid = linspace(0,6,1000);
nprobs = zeros(ni,length(xgrid));


for a = 1:ni
    % Each iteration, choose values for each min-max pair
    thismin = randn(size(mins)).*dmins + mins;
    thismax = randn(size(mins)).*dmaxs + maxs;
    % Everything between the highest min and the lowest max has equal
    % probability. 
    OK = ((xgrid <= min(thismax)) & (xgrid >= max(thismin))); % logical array
    this_pdf = double(OK);
    % Normalize to total probability = 1
    if any(OK)
        this_pdf = this_pdf./sum(this_pdf);
    end
    nprobs(a,:) = this_pdf;
end

% Summate all the iterations to get a summary distribution. 
total_prob = sum(nprobs)./ni;

% approximate normal distribution

weights = total_prob(total_prob > 0)./sum(total_prob(total_prob > 0));
values =  xgrid(total_prob>0);
mu = sum(values.*weights)./sum(weights);
sig = std(values, weights);


pdf_est = (1./(sig.*sqrt(2.*pi))).*exp(-0.5.*((xgrid - mu)./sig).^2);

% add Eaves distribution
% pd_eaves = makedist('Normal', 'mu', 3.6, 'sigma', 0.8);
% pdf_est_eaves = pdf(pd_eaves, xgrid);

subplot(2,1,2);
hold on 
stairs(xgrid,total_prob./sum(total_prob), 'k', 'LineWidth', 1.5);
plot(xgrid, pdf_est./sum(pdf_est), 'r', 'LineWidth', 1.5)
% plot(xgrid, pdf_est_eaves./sum(pdf_est_eaves), 'color', [0 0 1 0.4], 'LineWidth', 1.5)
txt = (['P_{10,sp} = ', sprintf('%0.1f', mu), ' ± ', sprintf('%0.1f', sig), ' at. g^{-1} yr^{-1}']);
text(4.2, 0.012, txt, 'FontSize', 14);

grid on;
box on;
set(gca,'xlim',[2 6], 'FontSize', 16, 'YTickLabel', []);
xlabel('P_{10,sp} (atoms g^{-1} yr^{-1})')
ylabel('Normalized Probability')
        
    
