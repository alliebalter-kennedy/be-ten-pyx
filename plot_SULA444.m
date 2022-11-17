% Plots results of leaching experiment for sample 444. 

samples.id = {'SULA-444-HF1', 'SULA-444-HF2', 'SULA-444-HF3', 'SULA-444-HF4', 'SULA-444-HF5'};

samples.N10 = [1.25E+07
1.06E+07
1.07E+07
1.04E+07
1.02E+07];

samples.dN10 = [2.96E+05
2.53E+05
3.38E+05
2.51E+05
2.46E+05];

samples.leaching_loss = [9;
25;
44;
52;
60];

avg = mean(samples.N10(2:end));
stdev = std(samples.N10(2:end));


labels = {'HF1', 'HF2', 'HF3', 'HF4', 'HF5'};
%% Plot

figure(1)
hold on

fill([0 70 70 0 0], [avg+stdev avg+stdev avg-stdev avg-stdev avg+stdev], 'k', 'EdgeColor','none', 'FaceAlpha',0.1)
plot([0 70], [avg avg], 'k--')
errorbar(samples.leaching_loss, samples.N10, samples.dN10, 'vertical', 'ks', 'MarkerFaceColor',[0.6 .6 0.6], 'MarkerSize',8)
text(samples.leaching_loss+1, samples.N10, labels,'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize',12)
set(gca, 'FontSize', 14)
xlabel('Sample Loss (%)', 'FontSize', 16)
ylabel('^{10}Be Concentration (atoms g^{-1})', 'FontSize', 16)
title('444 Leaching Experiment', 'FontSize', 16)
grid on

