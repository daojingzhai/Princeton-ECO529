function [pdf_stat,varargout] = KFE(X,MU,S,T,F0)
% KFE Solve one dimensional stationary and time-dependent KFE by finite
% difference method The process being studied is dX = MU(X)dt + S(X)dZ_t.

% REQUIRED INPUT:
% X: [X(1), X(2) ... X(N)]' is the state space (can be uneven grid)
% MU: a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% S: a volatility vector of lenght N, with S(1) = S(N) = 0
% OPTIONAL INPUT FOR TIME-DEPENDENT KFE: 
% T: time grid, M*1
% F0: initial distribution vector, N*1

% INPLEMENT
% 1. For stationary distribution, [pdf_stat] = KFE(X,MU,S);
% 2. For distribution diffusion, [pdf_stat,cdf] = KFE(X,MU,S,T,F0);

% NOTE: 1. Fokker–Planck operator (KFE) is the adjoint operator of Feynman–Kac 
% operator (KBE). We first build Feynman–Kac operator and then transpose it.
% 2. We use upwind scheme and implicit scheme for monotonicity and stability.

N = length(X);
dX = X(2:N) - X(1:N-1);
%% 1. Build Fokker-Planck operator
% approximate drift terms with an upwind scheme
% upper diagonal
AU = max(MU(1:N-1),0)./dX; 
% lower diagonal
AD = - min(MU(2:N),0)./dX;
% main diagonal 
A0 = zeros(N,1); A0(1:N-1) = A0(1:N-1) - AU; A0(2:N) = A0(2:N) - AD;
% matrix A
A = sparse(1:N,1:N,A0,N,N) + sparse(1:N-1,2:N,AU,N,N) + sparse(2:N,1:N-1,AD,N,N);

% approximate volatility terms
% sigma^2/(x_{n+1} - x_{n-1})
S0 = zeros(N,1); S0(2:N-1) = S(2:N-1).^2./(dX(1:N-2) + dX(2:N-1));
% upper diagonal
BU = S0(1:N-1)./dX;
% lower diagonal
BD = S0(2:N)./dX;
% main diagonal 
B0 = zeros(N,1); B0(1:N-1) = B0(1:N-1) - BU; B0(2:N) = B0(2:N) - BD;
% matrix B
B = sparse(1:N,1:N,B0,N,N) + sparse(1:N-1,2:N,BU,N,N) + sparse(2:N,1:N-1,BD,N,N);

% Fokker–Planck operator
FP = (A+B)'; 

%% 2. Find stationary distribution.
% MATLAB doesn't have build-in kernel solver for sparse matrix, for higher
% efficiency one can use online package like spnull, etc.
F_stat = null(full(FP)); 
cdf_stat = cumsum(F_stat(:,1)./sum(F_stat(:,1)));
pdf_stat = [0;(cdf_stat(2:end)-cdf_stat(1:end-1))./dX];

%% 3. Solve the time-dependent KFE.
if nargin == 5
    F = F0;
    DT = [0 T(2:end)-T(1:end-1)];
    pdf_diffusion = zeros(length(T),length(F));

    for i = 1:length(T)
        F = (speye(N,N) - DT(i)*FP)\F;
        pdf_diffusion(i,:) = F;
    end
    varargout{1} = pdf_diffusion;
end

end

