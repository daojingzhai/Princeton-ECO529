% This file computes the equilibrium in the setting of Brunnermeier and
% Sannikov, under assumptions that (1) experts can issue debt and equity,
% but must retain at least fraction chi_ of risk and (2) both experts and
% households have CRRA utility with relative risk aversion gamma
% the investment function is assumed to be of the form Phi(iota) = 
% log(kappa*iota + 1)/kappa, but this can be changed easily.
% 
% Platform: MATLAB R2019a
% script required: BasicFig.m, DistImpulseRespon.m, DistImpulseRespon_degen.m, MakeMovie.m
% function required: payoff_policy_growth.m, KFE.m, PercentileLine.m, FanChart.m

close all; clear; clc

%% Step 0: initialize parameters
a_e = .11; a_h = .03;       % production rate
rho_e = .06; rho_h = .05;   % discount rate
delta = .05; sigma = .10;   % decay rate/volatility
gamma = 2; phi = 10;        % CRRA param/investment function param
alpha = .5;                 % least equity fraction of experts, (0,1]

%% Step 1: Set Eta grid and terminal condition, initialize variables
% Eta Grid: uneven grid, with more points near 0 and 1
N = 1000;
zz = linspace(0.001,0.999,N)';
Eta = 3*zz.^2  - 2*zz.^3; 
dEta = Eta(2:N) - Eta(1:N-1);  

% variables
Q = ones(N,1);  	% price of capital q
Qp = zeros(N,1); 	% q'(eta)
SSQ = zeros(N,1);  	% sigma + sigma^q
Psi_e = zeros(N,1); % capital fraction of experts psi^e
Chi_e = zeros(N,1);	% risk fraction of experts chi^e
S = zeros(N,1);     % sigma^{eta}*eta
MU = zeros(N,1);  	% mu^{eta}*eta
MV_e = zeros(N,1); 	% mu^{v_e}
MV_h = zeros(N,1);	% mu^{v_h}
SV_e = zeros(N,1);	% sigma^{v_e}
SV_h = zeros(N,1);	% sigma^{v_h}
G = zeros(N,1); 	% payoff flow, zeros in chapter 3' setting

% terminal condition of value function
V_e = a_e^(-gamma)*Eta.^(1-gamma);         
V_h = a_e^(-gamma)*(1-Eta).^(1-gamma);

%% Step 2: iterative method
for t = 1:400
    %% Step 2(a): compute v' 
    GS_e = (Eta./V_e).^(1/gamma); 
    GS_h = ((1 - Eta)./V_h).^(1/gamma);
    GS  =  GS_e + GS_h;   % temporary variable for goods market clearing
    
    Vp_e = (V_e(2:N) - V_e(1:N-1))./dEta;   % \partial{v_e}
    Vp_h = (V_h(2:N) - V_h(1:N-1))./dEta;   % \partial{v_h}
    Vpl_e = Vp_e./V_e(1:N-1);               % \partial{v^e}/v^e 
    Vpl_h = Vp_h./V_h(1:N-1);               % \partial{v^h}/v^h 
    
    VVlp = Vpl_h - Vpl_e + 1./(Eta(1:N-1).*(1 - Eta(1:N-1)));
    
    %% Step 2(b): solve the equilibrium conditions via Newton's method
    % Step 2(b)0: Find q(0) in autarky economy (psi=0), use fzero function
    % to find root
    psi_e = 0; ssq = sigma; 
    q = fzero(@(x) a_h - (x-1)/phi - GS_h(1)*x^(1/gamma), [0 a_h*phi+1]);
    Q(1) = q; q_old = q;
    
    % Step 2(b)i: find q, psi, chi and sigma + sigma^q when psi<1  
    for n = 2:N 
        % errors given guess   
        ER = [log(q)/gamma + log(GS(n)) - log(a_e*psi_e + a_h*(1-psi_e) - (q-1)/phi);
              ssq*(q - (q - q_old)*(alpha*psi_e - Eta(n))/dEta(n-1)) - sigma*q;  
              a_e - a_h -  q*alpha*(alpha*psi_e - Eta(n))*ssq^2*VVlp(n-1)];  
        % matrix of derivatives of errors (could shorten it since q_old = q)
        QN = zeros(3,3); 
        
        QN(1,:) = [1/(q*gamma) + 1/((a_e - a_h)*psi_e + a_h - (q - 1)/phi)/phi, ...
          -(a_e - a_h)/((a_e - a_h)*psi_e + a_h - (q - 1)/phi), 0];   
      
        QN(2,:) = [ssq*(1 - (alpha*psi_e - Eta(n))/dEta(n-1)) - sigma, ...
          -ssq*(q-q_old)*alpha/dEta(n-1), q - (q-q_old)*(alpha*psi_e - Eta(n))/dEta(n-1)];  
      
        QN(3,:) = [- alpha*(alpha*psi_e - Eta(n))*ssq^2*VVlp(n-1), ...
          -q*alpha^2*ssq^2*VVlp(n-1), -2*q*alpha*(alpha*psi_e - Eta(n))*ssq*VVlp(n-1)];        
      
        % iterate based on Newton's method
        EN = [q; psi_e; ssq] - QN\ER;
        
        % if the boundary of the crisis regime has been reached. psi_e = 1 from now on 
        if EN(2) > 1
            break;
        end  
        
        % new guesses
        q = EN(1); psi_e = EN(2); ssq = EN(3);  
        
        % save results 
        Q(n) = EN(1); Psi_e(n) = EN(2); SSQ(n) = EN(3); 
        Qp(n) = (Q(n) - q_old)/dEta(n-1); q_old = EN(1);
        
    end
    
	% Step 2(b)ii: find q, psi, chi and sigma + sigma^q when psi=1
    n1 = n;
    for n = n1:N
        ER = log(q)/gamma + log(GS(n)) - log(a_e - (q - 1)/phi);        
        QN = 1/(q*gamma)  + 1/(a_e - (q - 1)/phi)/phi;
        EN = q - ER/QN;

        q = EN; Q(n) =  EN; Psi_e(n) = 1; Qp(n) = (EN - q)/dEta(n-1);
        SSQ(n) = 1./(1 - (max(alpha*Psi_e(n),Eta(n)) - Eta(n))*Qp(n)/Q(n))*sigma;
    end
    
    %% Step 2(c): compute drift and volatility of Eta and V
    Chi_e     = max(alpha*Psi_e,Eta);
    S         = (Chi_e - Eta).*SSQ; S(N) = 0;
    Iota      = (Q - 1)/phi; 
    Phi       = log(Q)/phi;
    A         = a_e*Psi_e + a_h*(1 - Psi_e) - Iota;
    
    CN_e      = GS_e.*Q.^(1/gamma - 1)./Eta; 
    CN_h      = GS_h.*Q.^(1/gamma - 1)./(1 - Eta);
    
    SV_e(2:N) = Vpl_e.*S(2:N);
    SV_h(2:N) = Vpl_h.*S(2:N); 
    
    VarSig_e  = -SV_e + S./Eta       + SSQ + (gamma - 1)*sigma;  
    VarSig_h  = -SV_h - S./(1 - Eta) + SSQ + (gamma - 1)*sigma;
   
    MU(2:N-1) = ((a_e - Iota(2:N-1))./Q(2:N-1) - CN_e(2:N-1)).*Eta(2:N-1) + ...
                S(2:N-1).*(VarSig_e(2:N-1) - SSQ(2:N-1)) + ...
                Eta(2:N-1).*SSQ(2:N-1).*(VarSig_e(2:N-1) - VarSig_h(2:N-1))*(1-alpha); 
            
    MV_e = rho_e - CN_e - (1 - gamma)*(Phi - delta - gamma*sigma^2/2 + SV_e*sigma);         
    MV_h = rho_h - CN_h - (1 - gamma)*(Phi - delta - gamma*sigma^2/2 + SV_h*sigma); 
    
    %% Step 2(d): update V_e and V_h 
    % lambda0 is dt*rho if dt is small, can be at most 1 (1 = policy iteration) 
    % it is more agressive to set lambda0 closer to 1, but code may not converge
    lambda0 = 0.8;   
    V_e = payoff_policy_growth(Eta, MV_e, MU, S, G, V_e, lambda0); 
    V_h = payoff_policy_growth(Eta, MV_h, MU, S, G, V_h, lambda0);
end

% Basic figure plotting
BasicFig;
% save data, and go to next step: 
save('Eta_S_MU','Eta','S','MU'); clear

%% Step 3: Distribution, density diffusion and animation
% simulate distributional impulse response, starting from stationary
% distribution
DistImpulseRespon;
% simulate distributional impulse response, starting from degenerated
% distribution
DistImpulseRespon_degen;
% simulate distributional impulse response (difference to unshocked path)
DistImpulseResponDiff

% make 2D/3D diffusion movie. note it may take >20 mins!!!
% MakeMovie;