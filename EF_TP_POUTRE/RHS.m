function B = RHS(f, L, N, M)
% RHS  Second membre B du systeme AU = B : B(i) = int f(x) w_i(x) dx
%   B(1) = B(N+2) = 0 (Dirichlet)

    h = L/(N+1);
    B = zeros(N+2, 1);

    w = @(x,xbar,h) ((x-(xbar-h))/h).*(x>xbar-h).*(x<xbar) + ((xbar+h-x)/h).*(x>xbar).*(x<xbar+h);

    for i = 2:N+1
        xbar  = (i-1)*h;
        wtemp = @(x) w(x,xbar,h);
        psi   = @(x) f.*wtemp(x);
        B(i)  = quadsimp(psi, xbar-h, xbar, M) + quadsimp(psi, xbar, xbar+h, M);
    end
end