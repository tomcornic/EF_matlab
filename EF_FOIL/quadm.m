function S = quadm(f,a,b,M)
% Integrale de f sur [a,b] par la methode du point milieu compose (M sous-intervalles)
% Utilisee avec M = 1 pour l'INTEGRATION REDUITE des termes de cisaillement
% (evite le "shear locking" quand S est grand devant EI/L^2).
dxt = (b-a)/M;
S = 0;
for j = 0:M-1
    xm = a + (j+0.5)*dxt;
    S = S + f(xm);
end
S = S*dxt;
end
