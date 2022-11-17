% This script looks at He-3 concentrations in Larsen/Farley paper and tries
% to figure out sigma and alpha. 

clear all; close all;

dd = [%4.35	0.15	3.9
4.54	0.28	6.1
4.27	0.14	8.6
3.6	0.22	13
3.13	0.16	16.8
2.29	0.01	24.1
1.74	0.07	29.2
1.61	0.07	34.4
%1.7	0.2	43.2
1.52	0.14	60.7
%4.34	0.14	3.9
2.28	0.09	34.4
%1.47	0.16	43.2
1.26	0.06	76.6
1	0.05	86.5
0.68	0.15	104.4
0.72	0.17	122.2
0.63	0.13	134.4
1	0.09	190.3
0.57	0.08	225.4
0.74	0.05	265.2
0.75	0.63	272.8
0.73	0.15	297.2
0.57	0.59	303.9];

d.sample_names = {%'CRB1IL'
'CRB2IL'
'CRB3IL'
'CRB4IL'
'CRB5IL'
'CRB6IL'
'CRB7IL'
'CRB8IL'
%'CRB9IL'
'CRB10IL'
%'CRB1PX'
'CRB8PX'
%'CRB9PX'
'CRB12PX'
'CRB13PX'
'CRB14PX'
'CRB15PX'
'CRB16PX'
'CRB18PX'
'CRB21PX'
'CRB22PX'
'CRB23PX'
'CRB24PX'
'CRB25PX'};

d.ilmenite = contains(d.sample_names,'IL');
d.pyroxene = ~d.ilmenite;

d.N3raw = dd(:,1).*1e6; d.dN3raw = dd(:,2).*1e6;
d.N3c(d.pyroxene) = d.N3raw(d.pyroxene);
d.N3c(d.ilmenite) = d.N3raw(d.ilmenite)./0.78;
d.dN3c(d.pyroxene) = d.dN3raw(d.pyroxene);
d.dN3c(d.ilmenite) = sqrt((d.dN3raw(d.ilmenite)./0.78).^2 + (0.02.*d.N3raw(d.ilmenite)./(0.78.^2)).^2);

d.N3c = d.N3c'; d.dN3c = d.dN3c';

d.zm = dd(:,3); % meters depth below basalt surface
d.zcm = d.zm.*100;
d.z = d.zcm.*2.71 + 180.*2; % total subsurface depth


%% Make a plot

figure; 
for a = 1:length(d.N3c)
    xx = d.N3c(a) + [-1 1].*d.dN3c(a);
    yy = [1 1].*d.z(a);
    plot(xx,yy,'k'); hold on;
    if d.ilmenite(a)
        plot(d.N3c(a),d.z(a),'ko','markerfacecolor',[1 0.5 0.5]);
    else
        plot(d.N3c(a),d.z(a),'ks','markerfacecolor',[0.5 0.5 1]);
    end
end

set(gca,'ydir','reverse'); grid on;


%% Do fitting


X0 = [8 1 1 1]; % X is sigma0, alpha, Nnuc, px/il conversion fx

opts = optimset('fmincon');
opts = optimset('display','iter');

% Note: px/il factor fixed at 0.78
[xopt,fval] = fmincon(@(X) objective_3_noerosion(X,d),X0,[],[],[],[],[0 0 0 0.78],[Inf Inf Inf 0.78],[],opts);

results = objective_3_noerosion(xopt,d,1);
    
figure; 
for a = 1:length(results.N3p)
    xx = results.N3obs(a) + [-1 1].*d.dN3c(a);
    yy = [1 1].*d.z(a);
    plot(xx,yy,'-','color',[0.5 0.5 0.5]); hold on; 
end
plot(results.N3obs,d.z,'o','color',[0.5 0.5 0.5],'markerfacecolor',[1 0.5 0.5]);

plot(results.N3p,d.z,'o','color','k','markerfacecolor','k');

set(gca,'ydir','reverse'); grid on;

set(gca,'xlim',[0 7e6])
xlabel('[^{3}He] (atoms g^{-1})');
ylabel('Mass depth (g cm^{-2})');
set(gca,'fontsize',12)

%% Now do with exponential fitting

X20 = [8 1 1 1];
[xopt2,fval2] = fmincon(@(X) objective_exponential(X,d),X0,[],[],[],[],[0 0 0 0.78],[Inf Inf Inf 0.78],[],opts);

L = xopt2(1).*1000;
N0 = xopt2(2).*1e6;
Nnuc = xopt2(3).*1e6

pz = 0:100:90000;
px = Nnuc + N0.*exp(-pz./L);
plot(px,pz,'b:');




