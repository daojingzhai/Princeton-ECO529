% This file simulates distributional impulse response (difference to
% unshocked paths). Suppose we have two systems staying in stationary
% distribution. We impose a shock (1% quantile of original Brownian shock,
% i.e. dZt = −2.32dt) for a period of dt = 1 to System 1, while don't
% impose it on System 2. After that, both systems experience a same
% sequence of Brownian shock (correlated Brownian). We simulate the density
% transition paths with Monte Carlo simulation and investigate their
% difference.
% Platform: MATLAB R2019a
% Data require: Eta_S_MU.mat 
% Funtion required: KFE.m, fanChartPath.m, PercentileLine.m

clc;clear
load('Eta_S_MU.mat')

Dt = 1;             % shock period
T0 = -10:0.5:-Dt;   % steady path grid
T1 = 0:0.1:Dt;      % shocked path grid
T2 = 0:1:200;       % transition path grid

E50 = [];   % 25/50/75 percentile line

%% Solve distribution
% solve stationary distribution
pdf_stat = KFE(Eta,MU,S);
sol0 = (pdf_stat*ones(1,length(T0)))';
cdf0 = cumsum((sol0(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol0,1))));
E50 = [E50 PercentileLine(T0, Eta, cdf0, 50)'];


%% shocked path
[~,sol1] = KFE(Eta,MU-2.32*S,0*S,T1,pdf_stat);
pdf_init = sol1(end,:)'; 
cdf1 = cumsum((sol1(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(sol1,1))));
E50 = [E50 PercentileLine(T1, Eta, cdf1, 50)'];

%% Monte-Carlo simulation

% number of path
n = 10000;
% number of time steps
m = 300;
% time step
h = 1;
% process
Eta1 = zeros(n,m);
Eta2 = zeros(n,m);

% starting point
eta1_0 = E50(1);
eta2_0 = E50(end);

% seed
rng(1)

for i = 1:n
    % starting point
    R = mvnrnd([eta1_0 eta2_0], [.000001 0; 0 .000001]);
    eta1 = R(1); eta2 = R(2);
    for j = 1:m
        rnum = normrnd(0,1);
        eta1 = eta1 + interp1(Eta,MU,eta1)*h + interp1(Eta,S,eta1)*sqrt(h)*rnum;
        eta2 = eta2 + interp1(Eta,MU,eta2)*h + interp1(Eta,S,eta2)*sqrt(h)*rnum;
        Eta1(i,j) = eta1;
        Eta2(i,j) = eta2;
    end
end

% save simulation results
% save('solution_DIR_Diff.mat');

% path differece
dEta = Eta2-Eta1;

%% fanchart
% 2D deta transition fan chart (colored by percentile, start from stationary
% distribution)
figure(8)
dE50 = E50-E50(1);
p0 = plot([T0 T1], dE50,'k','LineWidth',2);hold on
p1 = FanChartPath(1:size(dEta,2), dEta');hold on
xlabel('Time $t$','Interpreter','LaTex');
ylabel('$d \eta^e$','Interpreter','LaTex');
xline(0,'--b','Shock');xline(1,'--b');
xlim([-10 150]); ylim([-.025 .01]);
pbaspect([2 1 1]);
xticks([-10 0 50 100 150])
yticks([-.02 -.01 0 .01])
legend(p1,'Median State','Interpreter','LaTex','Location','southeast')


