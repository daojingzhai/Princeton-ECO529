function [sol] = KFE_pdepe(X,MU,S,tspan,intd)
% DynamicDist Solve time-dependent KFE to get diffusion path
% The process being studied is dX = MU(X)dt + S(X)dZ_t

% X = [X(1), X(2) ... X(N)]' is the state space (an increasing grid)
% MU is a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% S is a volatility vector of lenght N, with S(1) = S(N) = 0
% tspan is a vector specifying the points at which a solution is requested
% intd is the initial distribution

sol = pdepe(0,@pde1,@pdeic,@pdebc,X,tspan);

    function [c,f,s] = pde1(x,t,u,DuDx)
        N = length(X);
        mu = interp1(X,MU,x); 
        v = interp1(X,S,x);
        if x < X(2)
            dv = (S(2)-S(1))/(X(2)-X(1));
        else
            dv = interp1(X(2:N),(S(2:end)-S(1:end-1))./(X(2:end)-X(1:end-1)),x);
        end
        c = 1;
        f = -u*mu + v*dv*u + 0.5*v^2*DuDx;
        s = 0;
    end

    function u0 = pdeic(x) % pdeic: PDE initial condition
        u0 = interp1(X,intd,x);
    end

    function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t) % pdebc: PDE boundary condition
        pl = ul; ql = 0; pr = ur; qr = 0;
    end
end

