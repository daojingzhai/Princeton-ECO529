% solves chashless vs. monetary model in lecture 6
% required script: payoff_policy_growth.m
close all; clear; clc

%% parameters
a = 0.15;           % productivity
rho = 0.03;         % decay rate
sigma = 0.1;        % aggregate volatility
sigmaIdio = 0.5;    % indio volatility
phi = 2;            % investment function parameter 
chibar = 0.8;       % risk claims upperbound of intermediaries
varphi = 0.4;       % diversification ability of intermediaries for idio risk 
delta = 0.03;       % decay rate
tol = 1e-3;         % tolerance
muB = 0.01;         % monetary policy: constant
sigmaB = 0;         % monetary policy: constant
lambda = 0.5;       % lambda = dt/(1+dt), dt is the time step

%% Grid
etaLength = 300;
eta = (linspace(0.001,.95,etaLength))';  % wealth share of intermediaries
S = zeros(etaLength,1);           
G = zeros(etaLength,1);

%% Cashless Economy -- closed form solution
% real price of capital
qKCashless = (1+phi*a)/(1+phi*rho)*ones(etaLength,1);
% normalized price of outside money
qBCashless = zeros(etaLength,1);
% Investment
iotaCashless = (a-rho)/(1+phi*rho)*ones(etaLength,1);
% risk share of intermediaries
chiCashless = min(eta*(sigma^2+sigmaIdio^2)...
    ./(sigma^2 + ((1-eta)*varphi^2 + eta)*sigmaIdio^2), chibar);
% drift of etaI
muEtaCashless = (chiCashless-eta).^3*sigma^2./eta.^2./(1-eta) ...
    + (1-eta)*sigmaIdio^2 ...
    .*((chiCashless./eta).^2*varphi^2 - ((1-chiCashless)./(1-eta)).^2);
% volatility of etaI
sigmaEtaCashless = (chiCashless - eta)./eta * sigma;
% Risk free rate
varthetaCashless = zeros(etaLength,1);
PhiCashless = log(qKCashless)/phi;
riskFreeRateCashless = rho + PhiCashless - delta - chiCashless./eta*(sigma^2+sigmaIdio^2);

%% Monetary Economy -- 
% risk share of intermediaries
chiMonetary = min(eta./((1-eta)*varphi^2 + eta), chibar);
% initial vartheta: start from steady state
varthetaSteadyStateMonetary = 1 - sqrt(rho)/sigmaIdio/varphi;
varthetaInitMonetary = varthetaSteadyStateMonetary*ones(length(eta),1);
varthetaMonetary = varthetaInitMonetary;

for i=1:1500
  % 1. compute updated coefficients
  muEtaMonetary = (1-eta) .* (1-varthetaMonetary).^2 * sigmaIdio^2 ...
    .*((chiMonetary./eta).^2*varphi^2 - ((1-chiMonetary)./(1-eta)).^2);
  muVarthetaMonetary = rho + muB - (1-varthetaMonetary).^2 * sigmaIdio^2 ...
    .*(varphi^2*chiMonetary.^2./eta + (1-chiMonetary).^2./(1-eta));     % money valuation equation

  % 2. PDE time step, call payoff_policy_growth.m function
  MU = muEtaMonetary.*eta;
  newVarthetaMonetary = payoff_policy_growth(eta, muVarthetaMonetary, MU, S, G, varthetaMonetary, lambda); 

  % 3. check convergence
  absChangeVartheta = abs(newVarthetaMonetary-varthetaMonetary)/lambda*(1-lambda);
  relChangeVartheta = absChangeVartheta./(abs(newVarthetaMonetary)+abs(varthetaMonetary))*2;
  maxRelChange = max(relChangeVartheta);
  if maxRelChange < tol
      break;
  end
  
  % 4. update vartheta
  varthetaMonetary = newVarthetaMonetary;
end
% real price of capital
qKMonetary = (1-varthetaMonetary).*(1 + phi*a)./(1 - varthetaMonetary + phi*rho);
% normalized price of outside money
qBMonetary = varthetaMonetary.*(1 + phi*a)./(1 - varthetaMonetary + phi*rho);
% Investment
iotaMonetary = ((1-varthetaMonetary)*a-rho)./(1-varthetaMonetary+phi*rho);
% drift of etaI
muEtaMonetary = (1-eta) .* (1-varthetaMonetary).^2 * sigmaIdio^2 ...
    .*((chiMonetary./eta).^2*varphi^2 - ((1-chiMonetary)./(1-eta)).^2);
% volatility of etaI
sigmaEtaMonetary = zeros(etaLength);
% Find risk free rate
PhiMonetary = log(qKMonetary)/phi;
muVarthetaMonetary = rho + muB - (1-varthetaMonetary).^2 * sigmaIdio^2 ...
    .*(varphi^2*chiMonetary.^2./eta + (1-chiMonetary).^2./(1-eta));     
muP = (1./varthetaMonetary + 1./(1 - varthetaMonetary + phi*rho)).*muVarthetaMonetary.*varthetaMonetary;
riskFreeRateMonetary = PhiMonetary - delta + muP - sigma^2;

%% plots figures
figure
subplot(2,3,1)
hold on
plot(eta,qKCashless,'LineWidth',1,'Color','#0072BD');
plot(eta,qBCashless,'--','LineWidth',1,'Color','#0072BD');
plot(eta,qKMonetary,'LineWidth',1,'Color','#D95319');
plot(eta,qBMonetary,'--','LineWidth',1,'Color','#D95319');
ylabel('q^K (solid) / q^B (dashed)')
xlabel('\eta^I')
subplot(2,3,2)
hold on
plot(eta,chiCashless,'LineWidth',1,'Color','#0072BD');
plot(eta,chiMonetary,'LineWidth',1,'Color','#D95319');
ylabel('\chi^I')
xlabel('\eta^I')
subplot(2,3,3)
hold on
% plot(eta,riskFreeRateCashless,'LineWidth',1.5,'Color','#0072BD');
% plot(eta,riskFreeRateMonetary,'LineWidth',1.5,'Color','#D95319');
% ylabel('r^f')
plot(eta,varthetaCashless,'LineWidth',1,'Color','#0072BD');
plot(eta,varthetaMonetary,'LineWidth',1,'Color','#D95319');
ylabel('\vartheta')
xlabel('\eta^I')
subplot(2,3,4)
hold on
plot(eta,iotaCashless,'LineWidth',1,'Color','#0072BD');
plot(eta,iotaMonetary,'LineWidth',1,'Color','#D95319');
ylabel('\iota')
xlabel('\eta^I')
subplot(2,3,5)
hold on
plot(eta,muEtaCashless.*eta,'LineWidth',1,'Color','#0072BD');
plot(eta,muEtaMonetary.*eta,'LineWidth',1,'Color','#D95319');
% ylim([-.1 0.05])
ylabel('\eta^I \mu^{\eta^I}')
xlabel('\eta^I')
subplot(2,3,6)
hold on
plot(eta,sigmaEtaCashless.*eta,'LineWidth',1,'Color','#0072BD');
plot(eta,sigmaEtaMonetary.*eta,'LineWidth',1,'Color','#D95319');
ylabel('\eta^I \sigma^{\eta^I}')
xlabel('\eta^I')