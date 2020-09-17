% This program solves for the model in Lec 5: A One-Sector Monetary Model
% with Idiosyncratic Risk, and gives Fig 5.2.

% housekeeping 
clear; close all; clc

% parameters 
a = 0.1; rho = 0.05; phi = 10; gamma = 1; sigma = 0;
muM = 0.05; sigmaM = 0.05;

% equilibrium 
varrho = rho;
s_min = sqrt((rho+muM+(1-gamma)*sigma*sigmaM - sigmaM^2)/gamma);
sigmat = s_min:0.01:5;

vartheta = 1 - sqrt((varrho+muM+(1-gamma)*sigma*sigmaM - sigmaM^2)/gamma)./sigmat;
q_K = (1+phi*a)*(1 - vartheta)./(1 - vartheta + phi*varrho);
q_M = (1+phi*a)*vartheta./(1 - vartheta + phi*varrho);
iota = ((1 - vartheta)*a - varrho)./(1 - vartheta + phi*varrho);

figure
plot([0 s_min],[0 0],'Color','b','LineWidth',1.5); hold on
plot([s_min 5],[0 0],'--','Color','b'); hold on
plot([0 s_min],(1+phi*a)/(1+phi*rho)*[1 1],'Color','r','LineWidth',1.5); hold on
plot([s_min 5],(1+phi*a)/(1+phi*rho)*[1 1],'--','Color','r'); hold on
p1 = plot(sigmat,q_M,'Color','b','LineWidth',1.5); hold on
p2 = plot(sigmat,q_K,'Color','r','LineWidth',1.5); hold on
xlabel('$\tilde{\sigma}$','Interpreter','LaTex');
xlim([0 3]);
ylim([0 3.5]);
xticks([0:0.5:5])
yticks([0 1 (1+phi*a)/(1+phi*rho) 2 3])
set(gca, 'YTickLabel', {'0','1','$\displaystyle \frac{1+\phi a}{1+\phi\varrho}$','2','3'}, 'TickLabelInterpreter', 'latex');
legend([p1, p2], {'$q^M$','$q^K$'},'Location','northwest','Interpreter','LaTex');
daspect([1 2 1])