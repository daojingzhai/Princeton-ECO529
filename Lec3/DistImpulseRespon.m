% This file simulates distributional impulse response. 
% we start from stationary distribution, impose a shock (1% quantile of
% original Brownian shock, i.e. dZt = -2.32dt) for a period of Dt = 1, then
% investigate the density transition period by solving kolmogorov forward
% equation. We plot fan chart and 2D/3D density transition figure for the
% response. 
% Platform: MATLAB R2019a 
% Data require: Eta_S_MU.mat 
% funtion required: KFE.m, Fanchart.m, PercentileLine.m

clc;clear
load('Eta_S_MU.mat')

Dt = 1;             % shock period
T0 = -10:0.5:-Dt;   % steady path grid
T1 = 0:0.1:Dt;      % shocked path grid
T2 = 0:1:200;       % transition path grid

E25 = []; E50 = []; E75 = [];   % 25/50/75 percentile line
%% Solve distribution
% solve stationary distribution
pdf_stat = KFE(Eta,MU,S);
sol0 = (pdf_stat*ones(1,length(T0)))';
cdf0 = cumsum((sol0(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol0,1))));
E25 = [E25 PercentileLine(T0, Eta, cdf0, 25)'];
E50 = [E50 PercentileLine(T0, Eta, cdf0, 50)'];
E75 = [E75 PercentileLine(T0, Eta, cdf0, 75)'];

%% shocked path
% sol1 = KFE_pdepe(Eta,MU-2.32*S,0*S,T1,pdf_stat);      % use build in pdepe function
[~,sol1] = KFE(Eta,MU-2.32*S,0*S,T1,pdf_stat);
pdf_init = sol1(end,:)';
cdf1 = cumsum((sol1(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol1,1))));
E25 = [E25 PercentileLine(T1, Eta, cdf1, 25)'];
E50 = [E50 PercentileLine(T1, Eta, cdf1, 50)'];
E75 = [E75 PercentileLine(T1, Eta, cdf1, 75)'];

%% transition path
% sol2 = KFE_pdepe(Eta,MU,S,T2,pdf_init);       % use build in pdepe function
[~,sol2] = KFE(Eta,MU,S,T2,pdf_init);

cdf2 = cumsum((sol2(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol2,1))));
E25 = [E25 PercentileLine(T2, Eta, cdf2, 25)'];
E50 = [E50 PercentileLine(T2, Eta, cdf2, 50)'];
E75 = [E75 PercentileLine(T2, Eta, cdf2, 75)'];

save('solution_DIR.mat');
%% one-std shock plot
% 2D density Transition 
figure(3)
for i = [5:5:35 40:20:size(sol2,1)]
    plot(Eta,sol2(i,:),'Color',[1/3 1/3 1/3]); hold on
end
p1 = plot(Eta,pdf_init,'b','LineWidth',2); hold on
p2 = plot(Eta,pdf_stat,'r','LineWidth',2); hold on
legend([p1 p2],{'Shocked Distribution','Stationary Distribution'},'Interpreter','LaTex');
xlim([0 1]); ylim([0 18]);
set(gca,'XTick',0:0.1:1); set(gca,'YTick',0:3:18);
xlabel('$\eta^e$','Interpreter','LaTex');
ylabel('$f(\eta^e,t)$','Interpreter','LaTex');

% 3D density Transition
figure(4)
N = length(Eta);
p0 = surf(Eta,fliplr(T0*20),sol0,'EdgeColor','none','FaceColor',[1/3,1/3,0],'FaceAlpha',0.3); hold on
p1 = surf(Eta,fliplr(-T1*20),sol1,'EdgeColor','none','FaceColor','magenta','FaceAlpha',0.3); hold on
p2 = surf(Eta,T2,sol2,'EdgeColor','none','FaceColor',[1/3,1/3,1/3],'FaceAlpha',0.3); hold on
set(gca,'Ydir','reverse')
p3 = plot3(Eta,ones(N)*T2(1),  pdf_init,'LineWidth',2,'Color','b'); hold on
p4 = plot3(Eta,ones(N)*T2(end),pdf_stat,'LineWidth',2,'Color','r'); hold on
p5 = plot3(Eta,ones(N)*(-T1(end)*20),pdf_stat,'LineWidth',2,'Color','r'); hold on
legend([p0 p1 p2 p3(1) p4(1)],{'Stationary Period','Shock Period','Transitional Period','Shocked Distribution','Stationary Distribution'},...
    'Location','best','Interpreter','LaTex');
xlabel('$\eta^e$','Interpreter','LaTex');
ylabel('Time $t$','Interpreter','LaTex');
zlabel('Density $f(\eta^e,t)$','Interpreter','LaTex');
yticks([-50 -20 0 50 100 150 200])
yticklabels({'-5','-1','0','50','100','150','200'})
ylim([-50 200])

% 2D eta transition fan chart (colored based on density, start from
% stationary distribution)
figure(5)
N = length(Eta);
surf(Eta,T0,sol0,'EdgeColor','none'); hold on
surf(Eta,fliplr(-T1),sol1,'EdgeColor','none'); hold on
surf(Eta,T2,sol2,'EdgeColor','none'); hold on
plot3(E25,[T0 fliplr(-T1) T2],10000*ones(size(E25)),'--','LineWidth',2,'Color',[0.6350, 0.0780, 0.1840]); hold on
plot3(E75,[T0 fliplr(-T1) T2],10000*ones(size(E25)),'--','LineWidth',2,'Color',[0.6350, 0.0780, 0.1840]); hold on
plot3(E50,[T0 fliplr(-T1) T2],10000*ones(size(E25)),'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840]); hold on
xlabel('$\eta^e$','Interpreter','LaTex');
ylabel('Time $t$','Interpreter','LaTex');
zlabel('Density $f(\eta^e,t)$','Interpreter','LaTex');
yticks([-10 0 25 50 75 100])
yticklabels({'-10','0','25','50','75','100'})
ylim([-10 100]); xlim([0.35 0.5]); 
view(90,90); 
yline(-1,'--b');yline(0,'--b','Shock');
set(gca,'Xdir','reverse')
map = [linspace(1,1,9)' linspace(1,0,9)' linspace(1,0,9)'];
colormap(map);colorbar;

% 2D eta transition fan chart (colored by percentile, start from stationary
% distribution)
figure(6)
T0 = -10:0.5:-1; sol0 = (pdf_stat*ones(1,length(T0)))';
p0 = FanChart(T0, Eta, sol0, 'pdf');hold on
p1 = FanChart(fliplr(-T1), Eta, sol1, 'pdf');hold on
p2 = FanChart(T2, Eta, sol2, 'pdf');
xlabel('Time $t$','Interpreter','LaTex');
ylabel('$\eta^e$','Interpreter','LaTex');
xline(-1,'--b','Shock');xline(0,'--b');
xlim([-10 100]); ylim([.25 .5]);
pbaspect([2 1 1]);
xticks([-10 0 25 50 75 100])
yticks([.2 .25 .3 .35 .4 .45 .5])
legend(p0,'Median State','Interpreter','LaTex','Location','southeast')