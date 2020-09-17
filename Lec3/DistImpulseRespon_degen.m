% This file simulates distributional impulse response. we start from a
% degenerated distribution, impose a shock (1% quantile of original
% Brownian shock, i.e. dZt = −2.32dt) for a period of dt = 1, then
% investigate the density transition period by solving kolmogorov forward
% equation. We plot fan chart for the response.
% Platform: MATLAB R2019a 
% Data require: Eta_S_MU.mat
% Funtion required: KFE.m, Fanchart.m, PercentileLine.m

clc;clear
load('Eta_S_MU.mat')

dt = 1;             % shock period
T0 = -10:0.5:-dt;   % steady path grid
T1 = 0:0.1:dt;      % shocked path grid
T2 = 0:1:200;       % transition path grid

%% shocked path
% find median state of starting time
pdf_stat = KFE(Eta,MU,S);
sol0 = (pdf_stat*ones(1,length(T0)))';
cdf0 = cumsum((sol0(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol0,1))));
E50 = PercentileLine(T0, Eta, cdf0, 50);

% find state after shock by solve ode
odefun = @(t,y) interp1(Eta,MU-2.32*S,y);
[t,sol1] = ode45(odefun,T1,E50(1));

N = length(Eta);
pdf_init = normpdf(Eta,sol1(end),0.01);

%% transition path
% sol2 = KFE_pdepe(Eta,MU,S,T2,pdf_init);       % build in pdepe function
[~,sol2] = KFE(Eta,MU,S,T2,pdf_init);

%% one-std shock plot
figure(7)
T0 = -10:0.5:-1; 
p0 = plot(T0, E50,'k','LineWidth',2);hold on
p1 = plot(fliplr(-T1),sol1,'k','LineWidth',2);hold on
p2 = FanChart(T2, Eta, sol2, 'pdf');
xlabel('Time $t$','Interpreter','LaTex');
ylabel('$\eta^e$','Interpreter','LaTex');
xline(-1,'--b','Shock');xline(0,'--b');
xlim([-10 100]); ylim([.3 .5]);
pbaspect([2 1 1]);
xticks([-10 0 25 50 75 100])
yticks([.2 .25 .3 .35 .4 .45 .5])
legend(p0,'Median State','Interpreter','LaTex','Location','southeast')