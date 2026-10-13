
%quadrature num
clear all
close all
clc


I = exp(1)-1,%
F=@(x)(exp(x));a=0;b=1;

M = 50

for k=1:length(M)
    rg = quadrg(F,a,b,M(k));
    pm = quadpm(F,a,b,M(k));
    tpz = quadtpz(F,a,b,M(k));
    simp = quadsimp(F,a,b,M(k));
    Erg(k)=abs(I-rg);
    Epm(k) = abs(I - pm);
    Etpz(k) = abs(I - tpz);
    Esimp(k) = abs(I - simp);
end

%subplot(2,1,1);
%plot( M, Erg,M, Epm, M, Etpz, M, Esimp);
%legend('Rg')

%subplot(2,1,2);
%plot(log(M),log(Erg),log(M), log(Epm), log(M), log(Etpz),log(M),log(Esimp));
%legend('Rg')
n = 10;
A = LHS(1,1,n,M);

B = RHS(1,1,n,M)

U = A\B

x = linspace(0, 1, numel(U));
V = @(x) 1 - cosh(x - 1/2) / cosh(1/2);

plot(x, U, 'o-')
hold on
plot(x, V(x), '-')
hold off
legend('U', 'V')
