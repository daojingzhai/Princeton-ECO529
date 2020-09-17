function F = payoff_policy_growth(X, R, MU, S, G, V, lambda)
% X, R, MU, S, and G are column vectors of the same length
% X = [X(1), X(2) ... X(N)]' is the state space (an increasing grid)
% R is the discount rate minus growth
% MU is a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% S is a volatility vector of lenght N, with S(1) = S(N) = 0
% G is a payoff flow of length N
% We are solving the value function equation R(X)*F(t,X) = G(X) + MU(X)*F_x(t,X) + S(X)^2*F_xx(t,X)/2 + F_t(t,X) 
% To solve the equation over a small time interval dt with continuation value function V, set lambda = dt/(1 + dt) (or simply lambda = dt)
% To get the stationary solution of this equation, set lambda = 1: then V has no effect on the solution
% Given value F(t + dt, X) = V, the value of F(t, X) is found through an implicit scheme, by solving the matrix equation (R - MU*D - S^2*DD/2 + I/dt)*F(t) = G + F(t+dt)/dt, where D denotes the first derivative operator (in the upwind direction) and DD denotes the second derivative operator
% Multiplying both sides by lambda = dt/(1 + dt), we obtain the actual equation solved, ((R - MU*D - S^2*DD/2)*lambda + I*(1 - lambda))*F(t) = G*lambda + F(t+dt)*(1 - lambda) where the matrix multiplying F(t) is the matrix which is inverted

if or(MU(1) < 0, MU(end) > 0)
    disp('error: not true that MU(1) >= 0 and MU(N) <= 0');
end
if or(S(1) ~= 0, S(end) ~= 0)
    disp('error: not true that S(1) and S(N) both zero');
end

N = length(X); 
dX = X(2:N) - X(1:N-1); 

S0 = zeros(N,1); S0(2:N-1) = S(2:N-1).^2./(dX(1:N-2) + dX(2:N-1));
DU = zeros(N,1); DU(2:N) = - (max(MU(1:N-1),0) + S0(1:N-1))./dX*lambda; 
DD = zeros(N-1,1); DD = - (max(-MU(2:N),0) + S0(2:N))./dX*lambda; 

D0 = (1 - lambda)*ones(N,1) + lambda*R; 
D0(1:N-1) = D0(1:N-1) - DU(2:N); D0(2:N) = D0(2:N) - DD;
A = spdiags(D0,0,N,N) + spdiags(DU,1,N,N) + spdiags(DD(1:N-1),-1,N,N);
F = A\(G*lambda + V*(1 - lambda)); 
