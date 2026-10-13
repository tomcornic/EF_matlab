function A = LHS(c, L, N, M)
h = L/(N+1);
A = zeros(N+2, N+2);
A(1,1)     = 1;   % Dirichlet : u0 = 0
A(N+2,N+2) = 1;   % Dirichlet : u_{N+1} = 0

w  = @(x,xbar,h) ((x-(xbar-h))/h).*(x>xbar-h).*(x<xbar) + ((xbar+h-x)/h).*(x>xbar).*(x<xbar+h);
dw = @(x,xbar,h) (1/h)*(x>xbar-h).*(x<xbar) - (1/h)*(x>xbar).*(x<xbar+h);

for i = 2:N+1                       % on remplit une ligne par nœud intérieur
    xbar = (i-1)*h;                 % position du nœud i

    wtemp  = @(x) w(x,xbar,h);      % brique i
    dwtemp = @(x) dw(x,xbar,h);     % sa dérivée

    % --- diagonale A(i,i) : deux cellules ---
    psi    = @(x) dwtemp(x).^2 + c.*wtemp(x).^2;
    A(i,i) = quadsimp(psi, xbar-h, xbar, M) + quadsimp(psi, xbar, xbar+h, M);

    % --- sur-diagonale A(i,i+1) puis sous-diagonale par symétrie ---
    if i <= N
        wnext  = @(x) w(x,xbar+h,h);   % brique i+1
        dwnext = @(x) dw(x,xbar+h,h);
        psi2     = @(x) dwtemp(x).*dwnext(x) + c.*wtemp(x).*wnext(x);
        A(i,i+1) = quadsimp(psi2, xbar, xbar+h, M);
        A(i+1,i) = A(i,i+1);
    end
end
end