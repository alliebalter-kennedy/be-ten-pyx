% this script estimates the percent pyroxene and percent feldspar in whole 
% Ferrar dolerite as well as pyroxene separates from BGC and CRPG. We
% assume the measured elemental concentrations in LDEO pyroxenes represent
% pure pyroxene. Plagioclase endmembers are Anorthite (from Subramaniam
% 1956) and Albite (from Kracek and Neuvonen,  1952). 

% I can't remember the next steps! but i did do this before! i remember i
% ended up with a matrix for each sample with the % each of pyroxene,
% anorthite and albite. but I can't remember how I optimized to get there.
% did i use the optimization toolbox?

% i think i had a vector for each endmember that was:
% sample = [%Fe %Al %Ca %Na %Mg]
% pct_minerals = [%pyroxene %albite %anorthite]
% each sample was %pyroxene + %albite + %anorthite = 1
% 
% %Pyx*%Fe + %Alb*%Fe + %An*%Fe = %Fe(sample)

% fxn -- sum([pct_mins]'*sample)

% minimize -- sum((fxn - sample).^2))

% Aeq = [1 1 1];
% beq = 

% objective function 

% X(1) = % pyroxene
% X(2) = % albite
% X(3) = % anorthite

% compositions with oxides
% pyroxene_comp = [13.38 2.06 16.70 0.12 16.83];
% albite_comp = [0.02 19.65 0 11.07 0.04];
% anorthite_comp = [0.08 36.18 19.37 0.22 0];
% 
% wr_comp = [9.04 14.97 12.32 1.46 8.06];
% bgc_comp = [16.21 1.61 13.75 0.17 15.83];
% crpg_comp = [13.79 5.08 14.86 0.43 14.15];

% LDEO vs. BGC - LDEO has higher Ca, so may be augite while bgc may be
% pigeonite 

% compositions with element percents
pyroxene_comp = [8.51	1.09	11.94	0.09	10.15];
albite_comp = [0.02	10.4	0.00	8.21	0.02];
anorthite_comp = [0.06	19.15	13.85	0.16	0.00];
magnetite_comp = [70 0 0 0 0];

wr_comp = [6.6 7.9 8.8 1.1 4.9];
bgc_comp = [11.3 0.9 9.8 0.1 9.5];
crpg_comp = [9.6 2.7 10.6 0.3 8.5];

x0 = [0.33 0.33 0.29 0.04];
Aeq = [1 1 1 1];
beq = 1;

fun = @(x) est_plag_obj(x, bgc_comp, pyroxene_comp, albite_comp, anorthite_comp, magnetite_comp);

[optx fval] = fmincon(fun, x0, [], [], Aeq, beq,[0 0 0], [Inf Inf Inf]);


%% function for optimization

function out = est_plag_obj(x, sample, pyroxene_comp, albite_comp, anorthite_comp, magnetite_comp)

    comp_pred = [x(1) x(2) x(3) x(4)] * [pyroxene_comp; albite_comp; anorthite_comp; magnetite_comp];

    miss = sum((comp_pred - sample).^2);

    out = miss;

end
