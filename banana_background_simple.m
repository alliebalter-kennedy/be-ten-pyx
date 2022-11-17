% This script creates the background for a Be-10/He-3 banana plot. 

function [bananaPlot ax1] = banana_background_simple(HeBeRatio, Lsp, l10, rho, erosion, age);

% [figHandle productionRates.constants.P10q_St]
%% Set up production rates
    
    elevation = 0;
    
    p.constants = bedrock_constants();
    
    % Atmospheric pressure at site
    p.siteP = antatm(elevation); % site air pressure
    
    % Define production rate info
    p.SFsp = 1; % scaling factor
    p.Lsp = Lsp;
    
    p.P3sp = p.constants.P3p_St .* p.SFsp; 
    p.P10sp = (p.P3sp ./ HeBeRatio); 

%% line of constant exposure

erosion = erosion.*rho;
    
N10_tt = ((p.P10sp./l10) .* (1-exp(-l10.*age)))./p.P10sp;
N3_tt = (p.P3sp .* age)./p.P3sp;

R_tt = N10_tt./N3_tt;

%% line of steady erosion

N10_ee = (p.P10sp./(l10 + (erosion./Lsp)))./p.P10sp;
N3_ee = ((p.P3sp .* Lsp) ./ erosion)./p.P3sp ;

R_ee = (N10_ee./N3_ee);

    %% lines of constant age and exposure
 
    age = [0e6:1e4:30e6];
    erosion = [0:1e-6:500e-6].*rho;
    
    for a = 1:length(erosion)
         for b = 1:length(age)
            if erosion(a) == 0
                N10_te(a, b) = ((p.P10sp./l10) .* (1-exp(-l10.*age(b))))./p.P10sp;
                N3_te(a,b) = (p.P3sp .* age(b))./p.P3sp;
            elseif age(b) == max(age)
                N10_te(a, b) = (p.P10sp./(l10 + (erosion(a)./Lsp)))./p.P10sp;
                N3_te(a,b) = ((p.P3sp .* Lsp) ./ erosion(a))./p.P3sp;
            else
                N10_te(a, b) = ((p.P10sp ./ (l10 + (erosion(a)./Lsp))) .* (1-exp(-(l10 + (erosion(a)./Lsp)).*age(b))))./p.P10sp ;
                N3_te(a, b) = ((p.P3sp ./ (erosion(a)./Lsp)) .* (1-exp(-(erosion(a)./Lsp).*age(b))))./p.P3sp;
            end
         end
    end

    R_te = N10_te./N3_te;

    %% plot
    
    bananaPlot = figure;

    % plot lines of constant erosion
    for i = 1:length(erosion)  
        if mod(round(erosion(i)./rho, 6), 5e-6) == 0 %& erosion(i)./rho.*1e6 < 60
            plot(N3_te(i, :), R_te(i, :), '-', 'color', [0, 0, 1, 0.4], 'LineWidth', 0.5)

        else
        end
        hold on
    end

    % plot lines of constant age
    

    for i = 1:length(age)
         if mod(age(i), 1e6) ==0 & age(i) > 0
               plot(N3_te(:, i), R_te(:, i), '-', 'color', [0, 0, 0, 0.4], 'LineWidth', 0.5)

         else 
         end
        hold on
    end
    plot(N3_tt, R_tt, 'k', 'LineWidth', 1.5, 'color', [0 0 1 0.7])
    plot(N3_ee, R_ee, 'k', 'LineWidth', 1.5, 'color', [0 0 0])

    % plot labels

    for i = 1:length(erosion)  
        if mod(round(erosion(i)./rho, 6), 5e-6) == 0 & erosion(i)./rho.*1e6 < 30
            txt_ee_y = sum([(max(R_te(i, find(3e6 <= age & age<=4e6)), [], 2).*(1/2)) (min(R_te(i, find(3e6 <= age & age<=4e6)), [], 2).*(1/2))], 2);
            [val this_col] = min(abs((R_te(i, :) - txt_ee_y')));
            txt_ee_x = N3_te(i, this_col);
            txt = sprintf('%0i', round(erosion(i)./rho.*1e6, 0));
            text(txt_ee_x, txt_ee_y, txt, 'BackgroundColor', [1 1 1], 'Margin', 0.001, 'FontSize', 10, 'color', [0 0 1], 'HorizontalAlignment','center')
        else
        end
    end

    for i = 1:length(age)
         if mod(age(i), 2e6) ==0 & age(i) > 0 & age(i) <11e6
           txt_tt_x = sum([(max(N3_te(:, i), [], 1).*(2/3))' (min(N3_te(:, i), [], 1).*(1/3))'], 2);
           [val this_row] = min(abs((N3_te(:, i) - txt_tt_x)));
           txt_tt_y = R_te(this_row, i);

           txt = sprintf('%0d', age(i)./1e6);
           text(txt_tt_x, txt_tt_y, txt, 'BackgroundColor', [1 1 1], 'Margin', 0.001, 'FontSize', 10)
               
         else
         end
    end
    

    grid on
    box on
    
    

%     txt_tt = sprintf()
    
    xlabel('[^{3}He]^{*} (years)', 'FontSize', 16)
    ylabel('[^{10}Be]^{*}/[^{3}He]^{*}', 'FontSize', 16)
    set(gca, 'xlim', [0e6 8e6], 'FontSize', 16)
    

    ax1 = gca;
    
    % inset with smaller banana
%     ax2 = axes('Position',[0.53 0.55 0.35 0.35]);
% 
%     hold on
% 
%     plot(N3_tt, R_tt, 'k', 'LineWidth', 1.5, 'color', [0 0 1 0.7])
%     plot(N3_ee, R_ee, 'k', 'LineWidth', 1.5, 'color', [0.4 0.4 0.4])
%     xlabel('[^{3}He]^{*} (years)', 'FontSize', 12)
%     ylabel('[^{10}Be]^{*}/[^{3}He]^{*}', 'FontSize', 12)
%     set(ax2, 'xlim', [0 8e6], 'fontsize', 12) 
%     grid on
%     box on
    
end 